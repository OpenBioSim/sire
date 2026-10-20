/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
  *
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published
  *  by the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  *
  *  For full details of the license please see the LICENSE file
  *  that should have come with this distribution.
  *
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include "linearap.h"
#include "nmatrix.h"

#include "SireError/errors.h"

#include "tostring.h"

#include <limits>

#include <QDebug>

using namespace SireMaths;

static int getMinimumRow(const NMatrix &costs, int j)
{
    double minval = costs(0, j);
    int min_i = 0;

    for (int i = 1; i < costs.nRows(); ++i)
    {
        if (costs(i, j) < minval)
        {
            minval = costs(i, j);
            min_i = i;
        }
    }

    return min_i;
}

static bool isTiny(double val)
{
    const double tiny = 1.0e-10;

    return val > -tiny and val < tiny;
}

static void check_linear_assignment(const NMatrix &costs, const QVector<int> &rows_to_columns,
                                    const QVector<int> &columns_to_rows, const QVector<double> &u,
                                    const QVector<double> &v)
{
    int dim = costs.nRows();

    BOOST_ASSERT(int(costs.nColumns()) == dim);
    BOOST_ASSERT(rows_to_columns.count() == dim);
    BOOST_ASSERT(columns_to_rows.count() == dim);
    BOOST_ASSERT(u.count() == dim);
    BOOST_ASSERT(v.count() == dim);

    if (dim == 0)
        return;

    else if (dim == 1)
    {
        if (rows_to_columns[0] != 0 or columns_to_rows[0] != 0)
            throw SireError::program_bug(QObject::tr("How on earth did we get the single cost assignment wrong???"),
                                         CODELOC);

        return;
    }

    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            double reduced_cost = costs(i, j) - u[i] - v[j];

            if (isTiny(reduced_cost))
                reduced_cost = 0;

            if (reduced_cost < 0)
                throw SireError::program_bug(QObject::tr("Negative reduced cost %1 (%2,%3)\n%4\n%5\n%6\n%7")
                                                 .arg(reduced_cost)
                                                 .arg(i)
                                                 .arg(j)
                                                 .arg(Sire::toString(rows_to_columns))
                                                 .arg(Sire::toString(u), Sire::toString(v))
                                                 .arg(costs.toString()),
                                             CODELOC);
        }
    }

    for (int i = 0; i < dim; ++i)
    {
        int j = rows_to_columns[i];

        double reduced_cost = costs(i, j) - u[i] - v[j];

        if (not isTiny(reduced_cost))
            throw SireError::program_bug(QObject::tr("Non null reduced cost %1\n%2\n%3\n%4\n%5")
                                             .arg(reduced_cost)
                                             .arg(Sire::toString(rows_to_columns))
                                             .arg(Sire::toString(u), Sire::toString(v))
                                             .arg(costs.toString()),
                                         CODELOC);
    }

    QVector<int> matched(dim, 0);

    for (int i = 0; i < dim; ++i)
    {
        int j = rows_to_columns[i];

        if (matched[j])
        {
            throw SireError::program_bug(
                QObject::tr("Column is matched more than once! %1\n%2").arg(j).arg(Sire::toString(rows_to_columns)),
                CODELOC);
        }
        else
            matched[j] = 1;
    }

    for (int i = 0; i < dim; ++i)
    {
        int j = rows_to_columns[i];

        if (i != columns_to_rows[j])
            throw SireError::program_bug(QObject::tr("Error in row solution i != j (%1, %2, %3)\n%4\n%5")
                                             .arg(i)
                                             .arg(j)
                                             .arg(columns_to_rows[j])
                                             .arg(Sire::toString(rows_to_columns), Sire::toString(columns_to_rows)),
                                         CODELOC);
    }

    for (int j = 0; j < dim; ++j)
    {
        int i = columns_to_rows[j];

        if (j != rows_to_columns[i])
            throw SireError::program_bug(QObject::tr("Error in column solution j != i (%1, %2, %3)\n%4\n%5")
                                             .arg(j)
                                             .arg(i)
                                             .arg(rows_to_columns[i])
                                             .arg(Sire::toString(rows_to_columns), Sire::toString(columns_to_rows)),
                                         CODELOC);
    }
}

namespace SireMaths
{

    /** Solve the linear assignment problem where the costs are contained
        in the square matrix 'costs', returning the assignment of each row
        to each column

        For example, the costs could represent the distances from three
        taxis (A, B and C) to three people (R, S and T). The problem is
        to find the combination of taxis to people that minimises the total
        distance (and thus total time) for the taxis to get to the people.
        If the distance from taxi A to person R is D_{A,R}, then the
        cost matrix would look like this;

        costs = / D_{A,R}  D_{A,S}  D_{A,T} \
                | D_{B,R}  D_{B,S}  D_{B,T} |
                \ D_{C,R}  D_{C,S}  D_{C,T} /

        with taxis represented using rows, and people on columns.

        This function returns the index of each column assigned to each row,
        so if the solution was A => T, B => R and C => S then this function
        will return the vector;

        rows_to_columns = ( 2, 0, 1 )

        (i.e. taxi[0] => person[2], taxi[1] => person[0], taxi[2] => person[1])

        This function uses the algorithm presented in;

        "A shortest augmenting path algorithm for dense and sparse linear
         assignment problems", R. Jonker and A. Volgenant,
         Computing, 38, 325-340, 1987

        and is based on the freeware code provided by A. Jonker
        at <http://www.assignmentproblems.com/linearAP.htm>

        LAP.CPP

        Copyright notice

        Linear Assignment Problem - Freeware
        Codes LAPJV and LAPMOD: Jonker and Volgenant (direct download)

        <http://www.assignmentproblems.com/codes/LAP_dense.zip>
    */
    QVector<int> solve_linear_assignment(const NMatrix &costs, bool check_result)
    {
        const int dim = costs.nRows();

        if (dim == 0)
            return QVector<int>();

        if (dim != int(costs.nColumns()))
            throw SireError::invalid_arg(QObject::tr("Cannot solve the linear assignment problem as this code "
                                                     "only supports dense square matricies, and you have passed "
                                                     "in a rectangular matrix of dimension %1 by %2.")
                                             .arg(costs.nRows())
                                             .arg(costs.nColumns()),
                                         CODELOC);

        QVector<int> rows_to_columns(dim, 0);
        QVector<int> columns_to_rows(dim, -1);

        if (dim == 1)
        {
            // simple solution
            rows_to_columns[0] = 0;
            return rows_to_columns;
        }

        QVector<double> u(dim, 0); // dual variable: row reduction number
        QVector<double> v(dim, 0); // dual variable: column reduction number

        QVector<int> matches(dim, 0); // counts how many times a row can be matched
        QVector<int> free(dim, 0);    // list of unassigned rows

        // three procedures;
        //  (1) column reduction
        //  (2) reduction transfer from unassigned to assigned rows
        //  (3) augmenting reduction of unassigned rows

        // (1) column reduction
        //
        // Loop over the columns in reverse order and assign the row with the minimum cost
        //
        for (int j = dim - 1; j >= 0; --j)
        {
            int min_i = getMinimumRow(costs, j);
            v[j] = costs(min_i, j);

            matches[min_i] += 1;

            if (matches[min_i] == 1)
            {
                rows_to_columns[min_i] = j;
                columns_to_rows[j] = min_i;
            }
            else
            {
                // this row is already assigned - we can't assign the column
                columns_to_rows[j] = -1;
            }
        }

        // qDebug() << "\n(1) rows_to_columns" << rows_to_columns;
        // qDebug() << "(1) columns_to_rows" << columns_to_rows;
        // qDebug() << "(1) u" << u;
        // qDebug() << "(1) v" << v;
        // qDebug() << "free" << free;

        // (2) reduction transfer
        //
        // Aim is to further reduce the assigned rows but it doesn't
        // affect the reduction sum
        //
        int nfree = 0;

        for (int i = 0; i < dim; ++i)
        {
            if (matches[i] == 0)
            {
                free[nfree] = i;
                ++nfree;
            }
            else if (matches[i] == 1)
            {
                // transfer reduction from rows that are assigned once
                int j1 = rows_to_columns[i];

                double min = std::numeric_limits<double>::max();

                for (int j = 0; j < dim; ++j)
                {
                    if (j != j1)
                    {
                        double reduced_cost = costs(i, j) - v[j];
                        min = qMin(min, reduced_cost);
                    }
                }

                v[j1] -= min;
            }
        }

        // qDebug() << "\n(2) rows_to_columns" << rows_to_columns;
        // qDebug() << "(2) columns_to_rows" << columns_to_rows;
        // qDebug() << "(2) u" << u;
        // qDebug() << "(2) v" << v;
        // qDebug() << "free" << free;

        // (3) augmentating row reduction
        //
        // Aim is to find augmenting paths starting from unassigned rows
        // - we do this twice
        //
        for (int twice = 1; twice <= 2; ++twice)
        {
            // how many free rows are there?
            int numfree = nfree;
            nfree = 0;

            int k = 0;

            int old_k = -1;
            int old_numfree = -1;

            int nrepeat = 0;

            while (k < numfree)
            {
                if (k == old_k and numfree == old_numfree)
                {
                    ++nrepeat;

                    if (nrepeat > 10000)
                    {
                        qDebug() << "The linear assignment problem has not converged!!!";

                        throw SireError::program_bug(QObject::tr("The linear assignment problem is not converging!!!"),
                                                     CODELOC);
                    }
                }
                else
                {
                    nrepeat = 0;
                    old_k = k;
                    old_numfree = numfree;
                }

                int i = free[k];
                ++k;

                // find minimum and second minimum reduced cost over columns
                double first_min = costs(i, 0) - v[0];
                double second_min = std::numeric_limits<double>::max();

                int j1 = 0;
                int j2 = -1;

                for (int j = 1; j < dim; ++j)
                {
                    double reduced_cost = costs(i, j) - v[j];

                    if (reduced_cost < second_min)
                    {
                        if (reduced_cost >= first_min)
                        {
                            second_min = reduced_cost;
                            j2 = j;
                        }
                        else
                        {
                            second_min = first_min;
                            first_min = reduced_cost;
                            j2 = j1;
                            j1 = j;
                        }
                    }
                }

                int i0 = columns_to_rows[j1];

                if (first_min < second_min)
                    // change the reduction of the first minimum column
                    // to increase the minimum reduced cost in the row
                    // to the second minimum
                    v[j1] -= (second_min - first_min);

                else if (i0 >= 0)
                {
                    // swap columns j1 and j2 as j2 may be unassigned
                    j1 = j2;
                    i0 = columns_to_rows[j2];
                }

                // (re)assign i to j1, possibly de-assigning an i0
                rows_to_columns[i] = j1;
                columns_to_rows[j1] = i;

                if (i0 >= 0)
                {
                    // minimum column j1 assigned earlier
                    if (first_min < second_min)
                        // put in current k, and go back to that k.
                        // continue augmenting path i - j1 with i0.
                        free[--k] = i0;
                    else
                        // no further augmenting reduction possible.
                        // store i0 in list of free rows for next phase.
                        free[nfree++] = i0;
                }
            }
        }

        // augment solution for each free row
        QVector<int> collist(dim, 0); // list of columns to be scanned
        QVector<double> d(dim, 0);    // 'cost-distance' in augmenting path calculation
        QVector<int> pred(dim, 0);    // row-predecessor of column in
                                      // augmenting/alternating path

        // qDebug() << "\n(3)" << rows_to_columns;
        // qDebug() << "(3)" << columns_to_rows;
        // qDebug() << "(3) u" << u;
        // qDebug() << "(3) v" << v;

        // qDebug() << "free" << free;
        // qDebug() << "numfree ==" << nfree;

        for (int f = 0; f < nfree; ++f)
        {
            int freerow = free[f];

            // Dijkstra shortest path algorithm.
            // runs until unassigned column added to shortest path tree.
            for (int j = 0; j < dim; ++j)
            {
                d[j] = costs(freerow, j) - v[j];
                pred[j] = freerow;
                collist[j] = j; // init column list.
            }

            int low = 0; // columns in 0..low-1 are ready, now none.
            int up = 0;  // columns in low..up-1 are to be scanned for current minimum, now none.
                         // columns in up..dim-1 are to be considered later to find new minimum,
                         // at this stage the list simply contains all columns

            int last = -1;
            int endofpath = -1;

            bool unassignedfound = false;

            double min = std::numeric_limits<double>::max();

            do
            {
                if (up == low) // no more columns to be scanned for current minimum.
                {
                    last = low - 1;

                    // scan columns for up..dim-1 to find all indices for which new minimum occurs.
                    // store these indices between low..up-1 (increasing up).
                    min = d[collist[up++]];

                    for (int k = up; k < dim; ++k)
                    {
                        int j = collist[k];
                        double h = d[j];

                        if (h <= min)
                        {
                            if (h < min) // new minimum.
                            {
                                up = low; // restart list at index low.
                                min = h;
                            }

                            // new index with same minimum, put on undex up, and extend list.
                            collist[k] = collist[up];
                            collist[up++] = j;
                        }
                    }

                    // check if any of the minimum columns happens to be unassigned.
                    // if so, we have an augmenting path right away.
                    for (int k = low; k < up; ++k)
                    {
                        if (columns_to_rows[collist[k]] < 0)
                        {
                            endofpath = collist[k];
                            unassignedfound = true;
                            break;
                        }
                    }
                }

                if (not unassignedfound)
                {
                    // update 'distances' between freerow and all unscanned columns,
                    // via next scanned column.
                    int j1 = collist[low];
                    low++;

                    int i = columns_to_rows[j1];
                    double h = costs(i, j1) - v[j1] - min;

                    for (int k = up; k < dim; ++k)
                    {
                        int j = collist[k];
                        double v2 = costs(i, j) - v[j] - h;

                        if (v2 < d[j])
                        {
                            pred[j] = i;
                            if (v2 == min) // new column found at same minimum value
                            {
                                if (columns_to_rows[j] < 0)
                                {
                                    // if unassigned, shortest augmenting path is complete.
                                    endofpath = j;
                                    unassignedfound = true;
                                    break;
                                }
                                // else add to list to be scanned right away.
                                else
                                {
                                    collist[k] = collist[up];
                                    collist[up++] = j;
                                }
                            }

                            d[j] = v2;
                        }
                    }
                }
            } while (!unassignedfound);

            // update column prices.
            for (int k = 0; k <= last; ++k)
            {
                int j1 = collist[k];

                v[j1] += d[j1] - min;
            }

            // reset row and column assignments along the alternating path.
            int ii;

            do
            {
                ii = pred[endofpath];
                columns_to_rows[endofpath] = ii;

                int j1 = endofpath;
                endofpath = rows_to_columns[ii];
                rows_to_columns[ii] = j1;
            } while (ii != freerow);
        }

        if (check_result)
        {
            for (int i = 0; i < dim; ++i)
            {
                int j = rows_to_columns[i];
                u[i] = costs(i, j) - v[j];
            }

            // qDebug() << "\n(final) rows_to_columns" << rows_to_columns;
            // qDebug() << "(final) columns_to_rows" << columns_to_rows;
            // qDebug() << "(final) u" << u;
            // qDebug() << "(final) v" << v;

            // qDebug() << "\nTotal cost ==" << total_cost << "\n";

            // qDebug() << "CHECK2" << brute_force_linear_assignment(costs);

            check_linear_assignment(costs, rows_to_columns, columns_to_rows, u, v);
        }

        return rows_to_columns;
    }

    /** Return the total cost of the arrangement in 'rows_to_columns' using the
        costs in the linear assignment costs matrix in 'costs' */
    double calculate_total_cost(const NMatrix &costs, const QVector<int> &rows_to_columns)
    {
        int dim = costs.nRows();

        if (dim != int(costs.nColumns()))
            throw SireError::invalid_arg(QObject::tr("Cannot solve the linear assignment problem as this code "
                                                     "only supports dense square matricies, and you have passed "
                                                     "in a rectangular matrix of dimension %1 by %2.")
                                             .arg(costs.nRows())
                                             .arg(costs.nColumns()),
                                         CODELOC);

        if (rows_to_columns.count() != dim)
            throw SireError::invalid_arg(QObject::tr("Cannot calculate the total cost as the number of rows in wrong "
                                                     "(%1 vs. %2)")
                                             .arg(rows_to_columns.count())
                                             .arg(dim),
                                         CODELOC);

        double total_cost = 0;

        const int *rows_to_columns_array = rows_to_columns.constData();

        for (int i = 0; i < dim; ++i)
        {
            total_cost += costs(i, rows_to_columns_array[i]);
        }

        return total_cost;
    }

    void visit(const NMatrix &costs, QVector<int> &rows_to_columns, int k, QVector<int> &min_rows_to_columns,
               double &mincost)
    {
        static int level = -2;
        level = level + 1;
        rows_to_columns[k] = level;
        const int N = rows_to_columns.count();

        if (level == N - 1)
        {
            // evaluate the price of this combination
            double cost = calculate_total_cost(costs, rows_to_columns);
            // qDebug() << "TESTING:" << rows_to_columns << "cost ==" << cost << mincost;

            if (cost < mincost)
            {
                min_rows_to_columns = rows_to_columns;
                mincost = cost;
                // qDebug() << " *** NEW MINIMUM ***";
            }
        }
        else
        {
            for (int i = 0; i < N; ++i)
            {
                if (rows_to_columns[i] == -1)
                    visit(costs, rows_to_columns, i, min_rows_to_columns, mincost);
            }
        }

        level = level - 1;
        rows_to_columns[k] = -1;
    }

    /** Solve the linear assignment problem using a brute force algorithm
        (*very* poor scaling and very inefficient) */
    QVector<int> brute_force_linear_assignment(const NMatrix &costs)
    {
        const int dim = costs.nRows();

        if (dim == 0)
            return QVector<int>();

        if (dim != int(costs.nColumns()))
            throw SireError::invalid_arg(QObject::tr("Cannot solve the linear assignment problem as this code "
                                                     "only supports dense square matricies, and you have passed "
                                                     "in a rectangular matrix of dimension %1 by %2.")
                                             .arg(costs.nRows())
                                             .arg(costs.nColumns()),
                                         CODELOC);

        QVector<int> rows_to_columns(dim, -1);

        if (dim == 1)
        {
            // simple solution
            rows_to_columns[0] = 0;
            return rows_to_columns;
        }

        double mincost = std::numeric_limits<double>::max();

        QVector<int> min_rows_to_columns(dim, -1);

        visit(costs, rows_to_columns, 0, min_rows_to_columns, mincost);

        return min_rows_to_columns;
    }

    /** Return the minimum possible total cost from the linear assignment
        costs matrix 'costs' */
    double get_lowest_total_cost(const NMatrix &costs)
    {
        return calculate_total_cost(costs, solve_linear_assignment(costs));
    }

} // end of namespace SireMaths
