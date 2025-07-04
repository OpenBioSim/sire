/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include <QDateTime>
#include <QElapsedTimer>
#include <QFile>
#include <QHash>
#include <QRegularExpression>
#include <QSet>
#include <QTextStream>

#include <tuple>

#include "SireIO/amberformat.h"
#include "SireIO/amberprm.h"
#include "SireIO/amberrst7.h"

#include "SireMol/element.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/connectivity.h"
#include "SireMol/mgname.h"
#include "SireMol/selector.hpp"

#include "SireMol/atomcutting.h"
#include "SireMol/atomeditor.h"
#include "SireMol/atomidx.h"
#include "SireMol/cgatomidx.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/molidx.h"
#include "SireMol/reseditor.h"
#include "SireMol/residuecutting.h"
#include "SireMol/select.h"
#include "SireMol/trajectory.h"

#include "SireMol/amberparameters.h"

#include "SireCAS/trigfuncs.h"

#include "SireMM/amberparams.h"
#include "SireMM/atomljs.h"
#include "SireMM/cljnbpairs.h"
#include "SireMM/internalff.h"
#include "SireMM/ljparameter.h"
#include "SireMM/lj1264parameter.h"

#include "SireVol/cartesian.h"
#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireMaths/maths.h"
#include "SireUnits/units.h"

#include "SireBase/findexe.h"
#include "SireBase/tempdir.h"
#include "SireBase/progressbar.h"
#include "SireBase/propertylist.h"

#include "SireIO/errors.h"

#include "SireMove/flexibility.h"
#include "SireMove/internalmove.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/unittest.h"

#include "SireError/errors.h"
#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMM;
using namespace SireMol;
using namespace SireMove;
using namespace SireMaths;
using namespace SireCAS;
using namespace SireSystem;
using namespace SireVol;
using namespace SireUnits;
using namespace SireStream;
using namespace SireBase;

static const RegisterMetaType<AmberPrm> r_parm;

const RegisterParser<AmberPrm> register_amberparm;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const AmberPrm &parm)
{
    writeHeader(ds, r_parm, 3);

    SharedDataStream sds(ds);

    sds << parm.flag_to_line << parm.int_data << parm.float_data << parm.string_data << parm.ffield
        << parm.lj_exceptions << parm.warns << parm.comb_rules
        << static_cast<const MoleculeParser &>(parm);

    return ds;
}

/** Function called to rebuild the excluded atom lists */
void AmberPrm::rebuildExcludedAtoms()
{
    const int nmols = this->nMolecules();
    const int natoms = this->nAtoms();

    if (nmols <= 0 or natoms <= 0)
        return;

    // number of excluded atoms per atom
    const auto nnb_per_atom = this->intData("NUMBER_EXCLUDED_ATOMS");
    // list of all excluded atoms
    const auto inb = this->intData("EXCLUDED_ATOMS_LIST");

    if (nnb_per_atom.count() != natoms)
    {
        throw SireIO::parse_error(QObject::tr("The number of excluded atoms per atom array is not equal to the "
                                              "number of atoms! Should be %1, but is %2!")
                                      .arg(natoms)
                                      .arg(nnb_per_atom.count()),
                                  CODELOC);
    }

    // get the number of atoms in each molecule
    if (molnum_to_atomnums.count() != nmols + 1) // molnum_to_atomnums has a fake idx 0 value
    {
        throw SireIO::parse_error(QObject::tr("Disagreement of the number of molecules. %1 versus %2")
                                      .arg(molnum_to_atomnums.count() - 1)
                                      .arg(nmols),
                                  CODELOC);
    }

    const auto nnb_per_atom_data = nnb_per_atom.constData();
    const auto inb_data = inb.constData();

    excl_atoms = QVector<QVector<QVector<int>>>(nmols);
    auto excl_atoms_data = excl_atoms.data();
    int nnbidx = 0;

    // we need to loop through the excluded atoms in atomnum order, not molidx order.
    // This is because a single molecule can be spread over several non-contiguous blocks
    // of atom numbers!
    QVector<QPair<int, int>> atomnum_to_molnum_atomidx(natoms + 1);
    auto atomnum_to_molnum_atomidx_data = atomnum_to_molnum_atomidx.data();

    const auto &const_molnum_to_atomnums = molnum_to_atomnums;

    // get the molecule number and index of the atom in that molecule for every
    // atom identified by its atom number
    for (int molnum = 1; molnum <= nmols; ++molnum)
    {
        const auto &atomnums_in_mol = const_molnum_to_atomnums.at(molnum);
        const auto atomnums_in_mol_data = atomnums_in_mol.constData();

        for (int iat = 0; iat < atomnums_in_mol.count(); ++iat)
        {
            const int atomnum = atomnums_in_mol_data[iat];
            atomnum_to_molnum_atomidx_data[atomnum] = QPair<int, int>(molnum, iat);
        }

        excl_atoms_data[molnum - 1] = QVector<QVector<int>>(atomnums_in_mol.count());
    }

    // now loop over the atoms in sequence, looking up their excluded atoms (as these
    // are arranged in atomnum sequence. Extract the excluded atoms for each atom,
    // placing the result into excl_atoms, which is indexed by molnum then atomidx in
    // molecule
    for (int atomnum = 1; atomnum <= natoms; ++atomnum)
    {
        int molnum = atomnum_to_molnum_atomidx_data[atomnum].first;
        int atomidx = atomnum_to_molnum_atomidx_data[atomnum].second;

        const int nnb = nnb_per_atom_data[atomnum - 1]; // atomnum is 1-indexed

        QVector<int> ex(nnb);
        auto ex_data = ex.data();

        for (int k = 0; k < nnb; ++k)
        {
            if (nnbidx >= inb.count())
            {
                throw SireIO::parse_error(
                    QObject::tr("Disagreement of the number of excluded atom nb terms! %1 versus %2")
                        .arg(nnbidx)
                        .arg(inb.count()),
                    CODELOC);
            }

            ex_data[k] = inb_data[nnbidx];
            nnbidx += 1;
        }

        excl_atoms_data[molnum - 1][atomidx] = ex; // molnum is 1-indexed
    }
}

/** Internal function used to return the start index for the bonds of each
    molecule, and the number of bonds in each molecule */
QVector<QVector<int>> indexBonds(const QVector<qint64> &bonds, const QVector<int> &atom_to_mol, const int nmols)
{
    QVector<QVector<int>> molbonds(nmols);
    auto molbonds_data = molbonds.data();

    const auto bonds_data = bonds.constData();
    const auto atom_to_mol_data = atom_to_mol.constData();

    for (int i = 0; i < bonds.count(); i += 3)
    {
        // format is atom0-atom1-parameter
        int mol0 = atom_to_mol_data[bonds_data[i] / 3];     // divide by three as index is
        int mol1 = atom_to_mol_data[bonds_data[i + 1] / 3]; // into the coordinate array

        if (mol0 != mol1)
            throw SireIO::parse_error(QObject::tr("Something went wrong as there is a bond between two different "
                                                  "molecules! (%1, %2)")
                                          .arg(mol0)
                                          .arg(mol1),
                                      CODELOC);

        molbonds_data[mol0].append(i);
    }

    return molbonds;
}

/** Internal function used to return the start index for the angles of each
    molecule, and the number of angles in each molecule */
QVector<QVector<int>> indexAngles(const QVector<qint64> &angs, const QVector<int> &atom_to_mol, const int nmols)
{
    QVector<QVector<int>> molangs(nmols);
    auto molangs_data = molangs.data();

    const auto angs_data = angs.constData();
    const auto atom_to_mol_data = atom_to_mol.constData();

    for (int i = 0; i < angs.count(); i += 4)
    {
        // format is atom0-atom1-atom2-parameter
        int mol0 = atom_to_mol_data[angs_data[i] / 3];     // divide by three as index is
        int mol1 = atom_to_mol_data[angs_data[i + 1] / 3]; // into the coordinate array
        int mol2 = atom_to_mol_data[angs_data[i + 2] / 3];

        if (mol0 != mol1 or mol0 != mol2)
            throw SireIO::parse_error(
                QObject::tr("Something went wrong as there is a angle between more than one different "
                            "molecule! (%1, %2, %3)")
                    .arg(mol0)
                    .arg(mol1)
                    .arg(mol2),
                CODELOC);

        molangs_data[mol0].append(i);
    }

    return molangs;
}

/** Internal function used to return the start index for the dihedrals of each
    molecule, and the number of dihedrals in each molecule */
QVector<QVector<int>> indexDihedrals(const QVector<qint64> &dihs, const QVector<int> &atom_to_mol, const int nmols)
{
    QVector<QVector<int>> moldihs(nmols);
    auto moldihs_data = moldihs.data();

    const auto dihs_data = dihs.constData();
    const auto atom_to_mol_data = atom_to_mol.constData();

    for (int i = 0; i < dihs.count(); i += 5)
    {
        // format is atom0-atom1-atom2-atom3-parameter
        int mol0 = atom_to_mol_data[dihs_data[i] / 3];     // divide by three as index is
        int mol1 = atom_to_mol_data[dihs_data[i + 1] / 3]; // into the coordinate array
        int mol2 = atom_to_mol_data[std::abs(dihs_data[i + 2] / 3)];
        int mol3 = atom_to_mol_data[std::abs(dihs_data[i + 3] / 3)];

        if (mol0 != mol1 or mol0 != mol2 or mol0 != mol3)
            throw SireIO::parse_error(
                QObject::tr("Something went wrong as there is a dihedral between more than one different "
                            "molecule! (%1, %2, %3, %4)")
                    .arg(mol0)
                    .arg(mol1)
                    .arg(mol2)
                    .arg(mol3),
                CODELOC);

        moldihs_data[mol0].append(i);
    }

    return moldihs;
}

/** Internal function used to return the start index for the cmaps of each
    molecule, and the number of cmaps in each molecule */
QVector<QVector<int>> indexCMAPs(const QVector<qint64> &cmaps, const QVector<int> &atom_to_mol, const int nmols)
{
    QVector<QVector<int>> molcmaps(nmols);
    auto molcmaps_data = molcmaps.data();

    const auto cmaps_data = cmaps.constData();
    const auto atom_to_mol_data = atom_to_mol.constData();

    for (int i = 0; i < cmaps.count(); i += 6)
    {
        // format is atom0-atom1-atom2-atom3-atom4-parameter
        int mol0 = atom_to_mol_data[cmaps_data[i]];     // DO NOT DIVIDE BY THREE
        int mol1 = atom_to_mol_data[cmaps_data[i + 1]]; // THIS IS THE RAW ATOM INDEX
        int mol2 = atom_to_mol_data[cmaps_data[i + 2]];
        int mol3 = atom_to_mol_data[cmaps_data[i + 3]];
        int mol4 = atom_to_mol_data[cmaps_data[i + 4]];

        if (mol0 != mol1 or mol0 != mol2 or mol0 != mol3 or mol0 != mol4)
            throw SireIO::parse_error(
                QObject::tr("Something went wrong as there is a cmap between more than one different "
                            "molecule! (%1, %2, %3, %4, %5)")
                    .arg(mol0)
                    .arg(mol1)
                    .arg(mol2)
                    .arg(mol3)
                    .arg(mol4),
                CODELOC);

        molcmaps_data[mol0].append(i);
    }

    return molcmaps;
}

/** Function called to rebuild the Bond Angle and Dihedral indicies */
void AmberPrm::rebuildBADIndicies()
{
    const int nmols = this->nMolecules();
    const int natoms = this->nAtoms();

    if (nmols <= 0 or natoms <= 0)
        return;

    // get the lookup table to go from atom index to molecule index
    const auto atom_to_mol = this->getAtomIndexToMolIndex();

    // now index the connectivity - find the start index and number of bonds/angles/dihedrals
    // for each molecule
    if (usesParallel())
    {
        tbb::parallel_invoke(
            [&]()
            { bonds_inc_h = indexBonds(this->intData("BONDS_INC_HYDROGEN"), atom_to_mol, nmols); },
            [&]()
            { bonds_exc_h = indexBonds(this->intData("BONDS_WITHOUT_HYDROGEN"), atom_to_mol, nmols); },
            [&]()
            { angs_inc_h = indexAngles(this->intData("ANGLES_INC_HYDROGEN"), atom_to_mol, nmols); },
            [&]()
            { angs_exc_h = indexAngles(this->intData("ANGLES_WITHOUT_HYDROGEN"), atom_to_mol, nmols); },
            [&]()
            { dihs_inc_h = indexDihedrals(this->intData("DIHEDRALS_INC_HYDROGEN"), atom_to_mol, nmols); },
            [&]()
            { dihs_exc_h = indexDihedrals(this->intData("DIHEDRALS_WITHOUT_HYDROGEN"), atom_to_mol, nmols); },
            [&]()
            { cmap_idxs = indexCMAPs(this->intData("CMAP_INDEX") + this->intData("CHARMM_CMAP_INDEX"),
                                     atom_to_mol, nmols); });
    }
    else
    {
        bonds_inc_h = indexBonds(this->intData("BONDS_INC_HYDROGEN"), atom_to_mol, nmols);
        bonds_exc_h = indexBonds(this->intData("BONDS_WITHOUT_HYDROGEN"), atom_to_mol, nmols);
        angs_inc_h = indexAngles(this->intData("ANGLES_INC_HYDROGEN"), atom_to_mol, nmols);
        angs_exc_h = indexAngles(this->intData("ANGLES_WITHOUT_HYDROGEN"), atom_to_mol, nmols);
        dihs_inc_h = indexDihedrals(this->intData("DIHEDRALS_INC_HYDROGEN"), atom_to_mol, nmols);
        dihs_exc_h = indexDihedrals(this->intData("DIHEDRALS_WITHOUT_HYDROGEN"), atom_to_mol, nmols);
        cmap_idxs = indexCMAPs(this->intData("CMAP_INDEX") + this->intData("CHARMM_CMAP_INDEX"), atom_to_mol, nmols);
    }
}

/** Function called to rebuild all of the LJ parameters */
void AmberPrm::rebuildLJParameters()
{
    lj_data.clear();

    if (pointers.count() < 20)
        return;

    const auto pointers_data = pointers.constData();

    const int ntypes = pointers_data[1]; // number of distinct atom types

    if (ntypes <= 0)
        return;

    const int nphb = pointers_data[19]; // number of distinct 10-12 hydrogen bond pair types

    lj_data = QVector<LJParameter>(ntypes);
    auto lj_data_array = lj_data.data();

    lj_exceptions.clear();

    const auto acoeffs = float_data.value("LENNARD_JONES_ACOEF");
    const auto bcoeffs = float_data.value("LENNARD_JONES_BCOEF");
    const auto ccoeffs = float_data.value("LENNARD_JONES_CCOEF");

    const auto hbond_acoeffs = float_data.value("HBOND_ACOEF");
    const auto hbond_bcoeffs = float_data.value("HBOND_BCOEF");

    const auto nb_parm_index = int_data.value("NONBONDED_PARM_INDEX");

    if (acoeffs.count() != bcoeffs.count() or acoeffs.count() != (ntypes * (ntypes + 1)) / 2)
    {
        throw SireIO::parse_error(QObject::tr("Incorrect number of LJ coefficients for the number of specified "
                                              "atom types! Should be %1 for %2 types, but actually have "
                                              "%3 LJ A-coefficients, and %4 LJ B-coefficients")
                                      .arg((ntypes * (ntypes + 1)) / 2)
                                      .arg(ntypes)
                                      .arg(acoeffs.count())
                                      .arg(bcoeffs.count()),
                                  CODELOC);
    }

    if (nb_parm_index.count() != ntypes * ntypes)
    {
        throw SireIO::parse_error(QObject::tr("Incorrect number of non-bonded parameter indicies. There should "
                                              "be %1 indicies for %2 types, but actually have %3.")
                                      .arg(ntypes * ntypes)
                                      .arg(ntypes)
                                      .arg(nb_parm_index.count()),
                                  CODELOC);
    }

    if (hbond_acoeffs.count() != nphb or hbond_bcoeffs.count() != nphb)
    {
        throw SireIO::parse_error(QObject::tr("Incorrect number of HBond parameters. There should be "
                                              "%1 such parameters, but the number of HBond A coefficients is "
                                              "%2, and the number of B coefficients is %3.")
                                      .arg(nphb)
                                      .arg(hbond_acoeffs.count())
                                      .arg(hbond_bcoeffs.count()),
                                  CODELOC);
    }

    const auto acoeffs_data = acoeffs.constData();
    const auto bcoeffs_data = bcoeffs.constData();
    const auto ccoeffs_data = ccoeffs.constData();
    const auto hbond_acoeffs_data = hbond_acoeffs.constData();
    const auto hbond_bcoeffs_data = hbond_bcoeffs.constData();
    const auto nb_parm_index_data = nb_parm_index.constData();

    auto build_lj = [&](int i)
    {
        // amber stores the A and B coefficients as the product of all
        // possible combinations. We need to find the values from the
        //  LJ_i * LJ_i values
        int idx = nb_parm_index_data[ntypes * i + i];

        if (idx < 0)
        {
            // this is a 10-12 parameter
            throw SireError::unsupported(QObject::tr("Sire does not yet support Amber Parm files that "
                                                     "use 10-12 HBond parameters."),
                                         CODELOC);
        }
        else
        {
            double acoeff = acoeffs_data[idx - 1];
            double bcoeff = bcoeffs_data[idx - 1];

            double sigma = 0;
            double epsilon = 0;

            // numeric imprecision means that any parameter with acoeff less
            // than 1e-10 is really equal to 0
            if (acoeff > 1e-10)
            {
                // convert a_coeff & b_coeff into angstroms and kcal/mol-1
                sigma = std::pow(acoeff / bcoeff, 1 / 6.);
                epsilon = pow_2(bcoeff) / (4 * acoeff);
            }

            lj_data_array[i] = LJParameter(sigma * angstrom, epsilon * kcal_per_mol);
        }
    };

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, ntypes), [&](tbb::blocked_range<int> r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                build_lj(i);
            } });
    }
    else
    {
        for (int i = 0; i < ntypes; ++i)
        {
            build_lj(i);
        }
    }

    // The build_lj function above only considers diagonal elements of the
    // NONBONDED_PARM_INDEX matrix. Here we loop over the off-diagonal elements
    // to check for 10-12 parameters, which are currently unsupported.

    // The matrix is symmetric, so perform a triangular loop over off-diagonal
    // elements.
    for (int i = 0; i < ntypes; ++i)
    {
        for (int j = i + 1; j < ntypes; ++j)
        {
            int idx = nb_parm_index_data[ntypes * i + j];

            if (idx < 0)
            {
                auto a = hbond_acoeffs_data[1 - idx];
                auto b = hbond_bcoeffs_data[1 - idx];

                if ((a > 1e-6) and (b > 1e-6))
                {
                    // this is a 10-12 parameter
                    throw SireError::unsupported(QObject::tr("Sire does not yet support Amber Parm files that "
                                                             "use 10-12 HBond parameters."),
                                                 CODELOC);
                }
            }
        }
    }

    // While most LJ parameters use combining rules to form the i,j pairs,
    // Amber parm files support custom (exception) LJ parameters for specific
    // pairs of atoms. This is increasingly used by researchers to control
    // interactions between molecules, or to support the 12-6-4 potential
    // (which doesn't use combining rules?). Here, we loop over all pairs
    // of parameters to find those that don't use combining rules, and
    // so should be treated as exceptions.
    const bool has_c_coeffs = not ccoeffs.isEmpty();

    bool use_arithmetic_combining_rules = false;
    bool use_geometric_combining_rules = false;

    for (int i = 0; i < ntypes; ++i)
    {
        const auto &lj_i = lj_data_array[i];

        // include i==j pair as we also need to check for c-coeffs
        for (int j = i; j < ntypes; ++j)
        {
            const auto &lj_j = lj_data_array[j];

            int idx = nb_parm_index_data[ntypes * i + j];

            if (idx > 0)
            {
                auto a = acoeffs_data[idx - 1];
                auto b = bcoeffs_data[idx - 1];

                // convert A and B into sigma and epsilon
                double sigma = 0;
                double epsilon = 0;

                // numeric imprecision means that any parameter with acoeff less
                // than 1e-10 is really equal to 0
                if (a > 1e-10)
                {
                    // convert a_coeff & b_coeff into angstroms and kcal/mol-1
                    sigma = std::pow(a / b, 1 / 6.);
                    epsilon = pow_2(b) / (4 * a);
                }

                auto lj_ij = LJParameter(sigma * angstrom, epsilon * kcal_per_mol);

                auto expect = lj_i.combine(lj_j, LJParameter::ARITHMETIC);

                bool is_exception = false;

                if (std::abs(epsilon) < 1e-6 and std::abs(expect.epsilon().value()) < 1e-6)
                {
                    // this is a LJ pair that involves a ghost or dummy atom
                    // It should not impact exceptions or combining rules
                }
                else if (std::abs(lj_ij.epsilon().value() - expect.epsilon().value()) <= 1e-6)
                {
                    if (std::abs(lj_ij.sigma().value() - expect.sigma().value()) > 1e-6)
                    {
                        // not arithmetic combining rules - try geometric
                        expect = lj_i.combine(lj_j, LJParameter::GEOMETRIC);

                        is_exception = std::abs(lj_ij.sigma().value() - expect.sigma().value()) > 1e-6;

                        if (not is_exception)
                        {
                            if (has_c_coeffs)
                                use_geometric_combining_rules = ccoeffs_data[idx - 1] != 0;
                            else
                                use_geometric_combining_rules = true;
                        }
                    }
                    else
                    {
                        use_arithmetic_combining_rules = true;
                    }
                }
                else
                {
                    is_exception = true;
                }

                if (has_c_coeffs)
                {
                    // we have C coefficients, so we need to check if this is an exception
                    // or not
                    if (std::abs(ccoeffs_data[idx - 1]) > 1e-6)
                    {
                        // this is an exception
                        is_exception = true;
                    }
                }

                if (is_exception)
                {
                    // this is an exception parameter
                    double a = acoeffs_data[idx - 1];
                    double b = bcoeffs_data[idx - 1];
                    double c = 0;

                    if (has_c_coeffs)
                    {
                        c = ccoeffs_data[idx - 1];
                    }

                    // already in internal units (kcal mol-1, Angstroms)
                    auto lj_exception = LJException(LJ1264Parameter(a, b, c));

                    if (not lj_exceptions.contains(i))
                    {
                        lj_exceptions.insert(i, QList<LJException>());
                    }

                    if (not lj_exceptions.contains(j))
                    {
                        lj_exceptions.insert(j, QList<LJException>());
                    }

                    lj_exceptions[i].append(lj_exception);
                    lj_exceptions[j].append(lj_exception.getPair());
                }
            }
        }
    }

    if (use_arithmetic_combining_rules and use_geometric_combining_rules)
    {
        warns.append(QObject::tr("The LJ parameters in this Amber Parm file use both arithmetic and geometric "
                                 "combining rules. Sire will use arithmetic combining rules for all LJ "
                                 "parameters."));
        use_geometric_combining_rules = false;
    }

    if (use_geometric_combining_rules)
    {
        this->comb_rules = "geometric";
    }
    else
    {
        this->comb_rules = "arithmetic";
    }
}

/** This function is called to build all of the CMAP terms from the data
 *  read in from the PRM file
 */
void AmberPrm::rebuildCMAPTerms()
{
    cmap_data.clear();

    // some confusion as both CHARMM_X and X names are used
    const auto cmap_count = int_data.value("CMAP_COUNT") + int_data.value("CHARMM_CMAP_COUNT");

    if (cmap_count.size() < 2)
        return;

    // the first value is the CMAP_TERM_COUNT. This is not read - and we don't need
    // it because it is implied when we actually read in the CMAP terms
    // const int cmap_term_count = cmap_count[0];

    // the second is the CMAP_TYPE_COUNT
    const int cmap_type_count = cmap_count[1];

    const auto cmap_resolution = int_data.value("CMAP_RESOLUTION") + int_data.value("CHARMM_CMAP_RESOLUTION");

    // this should have CMAP_TERM_COUNT entries
    if (cmap_resolution.count() != cmap_type_count)
    {
        throw SireIO::parse_error(QObject::tr("The number of CMAP resolution terms is not equal to the number of "
                                              "CMAP terms! Should be %1, but is %2!")
                                      .arg(cmap_type_count)
                                      .arg(cmap_resolution.count()),
                                  CODELOC);
    }

    // there are CMAP_TYPE_COUNT data entries for the grids,
    // called CMAP_PARAMETER_01 to CMAP_PARAMETER_{CMAP_TYPE_COUNT}
    for (int i = 1; i <= cmap_type_count; ++i)
    {
        const auto resolution = cmap_resolution[i - 1];

        const auto cmap_parameter = float_data.value(QString("CMAP_PARAMETER_%1").arg(i, 2, 10, QChar('0'))) +
                                    float_data.value(QString("CHARMM_CMAP_PARAMETER_%1").arg(i, 2, 10, QChar('0')));

        // the number of values should equal the resolution squared
        if (cmap_parameter.count() != resolution * resolution)
        {
            throw SireIO::parse_error(QObject::tr("The number of CMAP parameters for type %1 is not equal to the "
                                                  "resolution squared! Should be %2, but is %3!")
                                          .arg(i)
                                          .arg(resolution * resolution)
                                          .arg(cmap_parameter.count()),
                                      CODELOC);
        }

        // the values are in column-major order
        cmap_data.insert(i, CMAPParameter(Array2D<double>::fromColumnMajorVector(
                                cmap_parameter, resolution, resolution)));
    }
}

/** This function finds all atoms that are bonded to the atom at index 'atom_idx'
(which is in molecule with index 'mol_idx', populating the hashe
    'atom_to_mol' (the molecule containing the passed atom). This uses the bonding information
    in 'bonded_atoms', which is the list of all atoms that are bonded to each atom */
static void findBondedAtoms(int atom_idx, int mol_idx, const QMultiHash<int, int> &bonded_atoms,
                            QHash<int, int> &atom_to_mol, QSet<int> &atoms_in_mol)
{
    for (auto bonded_atom : bonded_atoms.values(atom_idx))
    {
        if (not atoms_in_mol.contains(bonded_atom))
        {
            // we have not walked along this atom before
            atom_to_mol[bonded_atom] = mol_idx;
            atoms_in_mol.insert(bonded_atom);
            findBondedAtoms(bonded_atom, mol_idx, bonded_atoms, atom_to_mol, atoms_in_mol);
        }
    }
}

/** This function uses the bond information in 'bonds_inc_h' and 'bonds_exc_h'
    to divide the passed atoms into molecules. This returns an array of
    the number of atoms in each molecule (same format as ATOMS_PER_MOLECULE) */
static QVector<qint64> discoverMolecules(const QVector<qint64> &bonds_inc_h, const QVector<qint64> &bonds_exc_h,
                                         int natoms)
{
    // first, create a hash showing which atoms are bonded to each other

    // NOTE: the atom numbers in the following arrays that describe bonds
    // are coordinate array indexes for runtime speed. The true atom number
    // equals the absolute value of the number divided by three, plus one.
    //
    //%FORMAT(10I8)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
    //   IBH    : atom involved in bond "i", bond contains hydrogen
    //   JBH    : atom involved in bond "i", bond contains hydrogen
    //   ICBH   : index into parameter arrays RK and REQ
    QMultiHash<int, int> bonded_atoms;

    for (int j = 0; j < bonds_exc_h.count(); j = j + 3)
    {
        int atom0 = bonds_exc_h[j] / 3 + 1;
        int atom1 = bonds_exc_h[j + 1] / 3 + 1;
        bonded_atoms.insert(atom0, atom1);
        bonded_atoms.insert(atom1, atom0);
    }

    for (int j = 0; j < bonds_inc_h.count(); j = j + 3)
    {
        int atom0 = bonds_inc_h[j] / 3 + 1;
        int atom1 = bonds_inc_h[j + 1] / 3 + 1;
        bonded_atoms.insert(atom0, atom1);
        bonded_atoms.insert(atom1, atom0);
    }

    // Then recursively walk along each atom to find all the atoms that
    // are in the same molecule
    int nmols = 0;

    QHash<int, int> atom_to_mol;

    QList<qint64> atoms_per_mol;

    for (int i = 1; i <= natoms; i++)
    {
        if (not atom_to_mol.contains(i))
        {
            QSet<int> atoms_in_mol;

            nmols += 1;
            atom_to_mol[i] = nmols;
            atoms_in_mol.insert(i);

            // Recursive walk
            findBondedAtoms(i, nmols, bonded_atoms, atom_to_mol, atoms_in_mol);

            // this has now found all of the atoms in this molecule. Add the
            // number of atoms in the molecule to atoms_per_mol
            atoms_per_mol.append(atoms_in_mol.count());
        }
    }

    return atoms_per_mol.toVector();
}

/** Rebuild the arrays that show which atoms are in which molecules */
void AmberPrm::rebuildMolNumToAtomNums()
{
    const int natoms = this->nAtoms();

    const auto bonds_exc_h = this->intData("BONDS_WITHOUT_HYDROGEN");
    const auto bonds_inc_h = this->intData("BONDS_INC_HYDROGEN");

    const auto bonds_exc_h_data = bonds_exc_h.constData();
    const auto bonds_inc_h_data = bonds_inc_h.constData();

    // first, create a hash showing which atoms are bonded to each other

    // NOTE: the atom numbers in the following arrays that describe bonds
    // are coordinate array indexes for runtime speed. The true atom number
    // equals the absolute value of the number divided by three, plus one.
    //
    //%FORMAT(10I8)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
    //   IBH    : atom involved in bond "i", bond contains hydrogen
    //   JBH    : atom involved in bond "i", bond contains hydrogen
    //   ICBH   : index into parameter arrays RK and REQ
    QMultiHash<int, int> bonded_atoms;

    for (int j = 0; j < bonds_exc_h.count(); j = j + 3)
    {
        int atom0 = bonds_exc_h_data[j] / 3 + 1;
        int atom1 = bonds_exc_h_data[j + 1] / 3 + 1;
        bonded_atoms.insert(atom0, atom1);
        bonded_atoms.insert(atom1, atom0);
    }

    for (int j = 0; j < bonds_inc_h.count(); j = j + 3)
    {
        int atom0 = bonds_inc_h_data[j] / 3 + 1;
        int atom1 = bonds_inc_h_data[j + 1] / 3 + 1;
        bonded_atoms.insert(atom0, atom1);
        bonded_atoms.insert(atom1, atom0);
    }

    // Then recursively walk along each atom to find all the atoms that
    // are in the same molecule
    int nmols = 0;

    QHash<int, int> atom_to_mol;

    // remove the 0th index as Amber is 1-indexed
    QVector<QVector<int>> new_molnum_to_atomnums;

    new_molnum_to_atomnums.append(QVector<int>());

    QList<qint64> atoms_per_mol;

    // this is how many atoms the amber file says are in each molecule...
    const auto a_per_m = int_data.value("ATOMS_PER_MOLECULE");
    const auto a_per_m_data = a_per_m.constData();

    for (int i = 1; i <= natoms; i++)
    {
        if (not atom_to_mol.contains(i))
        {
            QSet<int> atoms_in_mol;

            nmols += 1;

            if (nmols > a_per_m.count())
            {
                throw SireIO::parse_error(
                    QObject::tr("The files appears to contain more molecules than expected. The file "
                                "should only contain %1 molecules, but more than this have been "
                                "found based on the bonding of the molecules.")
                        .arg(nmols),
                    CODELOC);
            }

            // the number of atoms we should expect in this molecule...
            int expected_nats = a_per_m_data[nmols - 1];

            // this is the first atom in the new molecule
            atom_to_mol[i] = nmols;
            atoms_in_mol.insert(i);

            // Recursive walk to find all of the other atoms
            findBondedAtoms(i, nmols, bonded_atoms, atom_to_mol, atoms_in_mol);

            // this has now found all of the atoms in this molecule. Check we have
            // all of the molecules expected by the amber file...
            int next_atom = i + 1;

            while (atoms_in_mol.count() < expected_nats)
            {
                // find the next atom which is not yet in a molecule...
                for (; next_atom <= natoms; ++next_atom)
                {
                    if (not atom_to_mol.contains(next_atom))
                    {
                        break;
                    }
                }

                // add this, and all of its bonded atoms
                atom_to_mol[next_atom] = nmols;
                atoms_in_mol.insert(next_atom);
                findBondedAtoms(next_atom, nmols, bonded_atoms, atom_to_mol, atoms_in_mol);

                next_atom += 1;
                if (next_atom > natoms)
                    break;
            }

            if (atoms_in_mol.count() != expected_nats)
            {
                throw SireIO::parse_error(
                    QObject::tr("Disagreement over the number of atoms in molecule %1. Looking at bonding "
                                "implies the number of atoms is %2, while the file itself claims the "
                                "number is %3.")
                        .arg(nmols)
                        .arg(atoms_in_mol.count())
                        .arg(expected_nats),
                    CODELOC);
            }

            // add the number of atoms in the molecule to atoms_per_mol
            atoms_per_mol.append(atoms_in_mol.count());

            auto atms = atoms_in_mol.values();
            std::sort(atms.begin(), atms.end());

            new_molnum_to_atomnums.append(atms.toVector());
        }
    }

    molnum_to_atomnums = new_molnum_to_atomnums;

    // now have all of the atomidxs (1-indexed) for molecule i (1-indexed)
    // in the array molidx_to_atomidxs. Guaranteed to be sorted into AtomNum order
    // in molnum_to_atomnums
}

/** Function called after loading the AmberPrm from a binary stream
    to populate all of the calculated member data */
void AmberPrm::rebuildAfterReload()
{
    pointers = this->intData("POINTERS");

    if (pointers.isEmpty())
    {
        this->operator=(AmberPrm());
    }
    else if (pointers.count() < 31)
    {
        throw SireIO::parse_error(
            QObject::tr("There was no, or an insufficient 'POINTERS' section in the file! (%1)").arg(pointers.count()),
            CODELOC);
    }

    // first need to create the map of all of the atom numbers in each molecule
    this->rebuildMolNumToAtomNums();

    if (usesParallel())
    {
        tbb::parallel_invoke(
            // now we have to build the LJ parameters (Amber stores them weirdly!)
            [&]()
            { this->rebuildLJParameters(); },
            // now we have to build the lookup indicies for the bonds, angles and dihedrals
            [&]()
            { this->rebuildBADIndicies(); },
            // now we have to build the excluded atom lists
            [&]()
            { this->rebuildExcludedAtoms(); },
            // now we have to build the CMAP terms
            [&]()
            {
                this->rebuildCMAPTerms();
            });
    }
    else
    {
        this->rebuildLJParameters();
        this->rebuildBADIndicies();
        this->rebuildExcludedAtoms();
        this->rebuildCMAPTerms();
    }
}

/** Read from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AmberPrm &parm)
{
    VersionID v = readHeader(ds, r_parm);

    if (v == 1 or v == 2 or v == 3)
    {
        SharedDataStream sds(ds);

        parm = AmberPrm();

        sds >> parm.flag_to_line >> parm.int_data >> parm.float_data >> parm.string_data;

        if (v == 2)
            sds >> parm.ffield;
        else
            parm.ffield = MMDetail();

        parm.lj_exceptions.clear();
        parm.comb_rules = "arithmetic";

        if (v == 3)
            sds >> parm.lj_exceptions >> parm.warns >> parm.comb_rules;
        else
            parm.lj_exceptions.clear();

        sds >> static_cast<MoleculeParser &>(parm);

        parm.rebuildAfterReload();
    }
    else
        throw version_error(v, "1,2,3", r_parm, CODELOC);

    return ds;
}

/** Constructor */
AmberPrm::AmberPrm() : ConcreteProperty<AmberPrm, MoleculeParser>()
{
    for (int i = 0; i < 32; ++i)
    {
        pointers.append(0);
    }
}

AmberPrm::FLAG_TYPE flagType(const QStringList &lines, const QPair<qint64, qint64> &index)
{
    AmberFormat f(lines[index.first + 1]);
    return f.flagType();
}

/** Return the flag type for the data associated with the passed flag.
    This returns UNKNOWN if this is not known */
AmberPrm::FLAG_TYPE AmberPrm::flagType(const QString &flag) const
{
    if (flag_to_line.contains(flag))
    {
        auto index = flag_to_line.value(flag);

        if (index.first > 0)
            return AmberFormat(lines()[index.first - 1]).flagType();
    }

    return AmberPrm::UNKNOWN;
}

/** Return the integer data for the passed flag. This returns an empty
    list if there is no data associated with this flag. This raises
    an invalid_cast error if data exists, but it is the wrong type */
QVector<qint64> AmberPrm::intData(const QString &flag) const
{
    auto it = int_data.constFind(flag);

    if (it != int_data.constEnd())
    {
        return it.value();
    }

    if (flag_to_line.contains(flag))
    {
        if (float_data.contains(flag))
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the float data for flag '%1' to integer data!").arg(flag), CODELOC);
        else if (string_data.contains(flag))
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the string data for flag '%1' to integer data!").arg(flag), CODELOC);
        else
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the data for flag '%1' to integer data!").arg(flag), CODELOC);
    }

    return QVector<qint64>();
}

/** Return the float data for the passed flag. This returns an empty
    list if there is no data associated with this flag. This raises
    an invalid_cast error if data exists, but it is the wrong type */
QVector<double> AmberPrm::floatData(const QString &flag) const
{
    auto it = float_data.constFind(flag);

    if (it != float_data.constEnd())
    {
        return it.value();
    }

    if (flag_to_line.contains(flag))
    {
        if (int_data.contains(flag))
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the integer data for flag '%1' to float data!").arg(flag), CODELOC);
        else if (string_data.contains(flag))
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the string data for flag '%1' to float data!").arg(flag), CODELOC);
        else
            throw SireError::invalid_cast(QObject::tr("Cannot convert the data for flag '%1' to float data!").arg(flag),
                                          CODELOC);
    }

    return QVector<double>();
}

/** Return the string data for the passed flag. This returns an empty
    list if there is no data associated with this flag. This raises
    an invalid_cast error if data exists, but it is the wrong type */
QVector<QString> AmberPrm::stringData(const QString &flag) const
{
    auto it = string_data.constFind(flag);

    if (it != string_data.constEnd())
    {
        return it.value();
    }

    if (flag_to_line.contains(flag))
    {
        if (float_data.contains(flag))
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the float data for flag '%1' to string data!").arg(flag), CODELOC);
        else if (int_data.contains(flag))
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the integer data for flag '%1' to string data!").arg(flag), CODELOC);
        else
            throw SireError::invalid_cast(
                QObject::tr("Cannot convert the data for flag '%1' to string data!").arg(flag), CODELOC);
    }

    return QVector<QString>();
}

/** Return the forcefield for the molecules in this file */
MMDetail AmberPrm::forcefield() const
{
    return ffield;
}

/** Process all of the flags */
double AmberPrm::processAllFlags()
{
    QMutex int_mutex, float_mutex, string_mutex;

    double score = 0;

    const auto const_flag_to_line = flag_to_line;

    const QVector<QString> flags = flag_to_line.keys().toVector();
    const auto flags_data = flags.constData();

    QStringList global_errors;

    const auto &l = lines();
    const auto lines_data = l.constData();

    auto process_flag = [&](int i, QStringList &errors, double &scr)
    {
        const QString &flag = flags_data[i];
        const QPair<qint64, qint64> index = const_flag_to_line.value(flag);

        // the format for the data is on the preceeding line
        const AmberFormat format(lines_data[index.first - 1]);

        switch (format.flagType())
        {
        case INTEGER:
        {
            QVector<qint64> data = readIntData(l, format, index, &scr, &errors);
            QMutexLocker lkr(&int_mutex);
            int_data.insert(flag, data);
            break;
        }
        case FLOAT:
        {
            QVector<double> data = readFloatData(l, format, index, &scr, &errors);
            QMutexLocker lkr(&float_mutex);
            float_data.insert(flag, data);
            break;
        }
        case FFLOAT:
        {
            QVector<double> data = readFloatData(l, format, index, &scr, &errors);
            QMutexLocker lkr(&float_mutex);
            float_data.insert(flag, data);
            break;
        }
        case STRING:
        {
            QVector<QString> data = readStringData(l, format, index, &scr, &errors);
            QMutexLocker lkr(&string_mutex);
            string_data.insert(flag, data);
            break;
        }
        default:
            break;
        }
    };

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, flags.count()), [&](tbb::blocked_range<int> r)
                          {
            QStringList local_errors;
            double local_score = 0;

            for (int i = r.begin(); i < r.end(); ++i)
            {
                process_flag(i, local_errors, local_score);
            }

            QMutexLocker lkr(&int_mutex);
            score += local_score;

            if (not local_errors.isEmpty())
            {
                global_errors += local_errors;
            } });
    }
    else
    {
        for (int i = 0; i < flags.count(); ++i)
        {
            process_flag(i, global_errors, score);
        }
    }

    if (not global_errors.isEmpty())
    {
        throw SireIO::parse_error(QObject::tr("There were errors parsing the file:%1").arg(global_errors.join("\n")),
                                  CODELOC);
    }

    // now have to get the general pointers array
    pointers = this->intData("POINTERS");

    if (pointers.isEmpty())
    {
        this->operator=(AmberPrm());
        return 0;
    }
    else if (pointers.count() < 31)
    {
        throw SireIO::parse_error(
            QObject::tr("There was no, or an insufficient 'POINTERS' section in the file! (%1)").arg(pointers.count()),
            CODELOC);
    }

    // now we have to work out how the atoms divide into molecules
    const auto atoms_per_mol = this->intData("ATOMS_PER_MOLECULE");

    if (atoms_per_mol.count() == 0)
    {
        auto atoms_per_mol =
            discoverMolecules(int_data.value("BONDS_INC_HYDROGEN"), int_data.value("BONDS_WITHOUT_HYDROGEN"), nAtoms());

        if (atoms_per_mol.count() == 0)
        {
            atoms_per_mol = QVector<qint64>(1, nAtoms());
        }

        int_data.insert("ATOMS_PER_MOLECULE", atoms_per_mol);
    }

    // some older prm7 files don't provide any information for ATOMIC_NUMBER. We need
    // to add dummy information here
    if (this->intData("ATOMIC_NUMBER").count() == 0)
    {
        auto atomic_number = QVector<qint64>(nAtoms(), 0);
        int_data.insert("ATOMIC_NUMBER", atomic_number);
    }

    // build everything that can be derived from this information
    this->rebuildAfterReload();

    return score;
}

/** Parse the text data from this file to create this object. This
    clears any pre-existing data */
void AmberPrm::parse(const PropertyMap &map)
{
    double score = 0;

    // as we are reading, look out for any FLAGs, so that we
    // can record their locations
    QString last_flag = QString();

    // make sure that there are no comment lines
    this->removeCommentLines("%COMMENT");

    const int nlines = lines().count();
    const QString *lines_array = lines().constData();

    for (int i = 0; i < nlines; ++i)
    {
        const QString &line = lines_array[i];

        if (line.length() > 0 and line[0] == '%')
        {
            // this is a control line
            if (line.startsWith("%FLAG"))
            {
                // We currently don't support CHAMBER files. Raise an exception
                // whenever a FLAG line contains CHARMM.
                if (line.contains("CHARMM"))
                {
                    throw SireError::file_error(QObject::tr("We currently don't support CHAMBER files."), CODELOC);
                }

                // this is a new flag - close any open old flag
                if (not last_flag.isNull())
                {
                    if (flag_to_line.contains(last_flag))
                    {
                        flag_to_line[last_flag].second = i - flag_to_line[last_flag].first;
                    }

                    last_flag = QString();
                }

                // find the new flag
                QStringList words = line.split(" ", Qt::SkipEmptyParts);

                QString flag = words[1];

                if (flag_to_line.contains(flag))
                    throw SireError::file_error(QObject::tr("The file does not look like a valid Amber Parm7 file, "
                                                            "as the FLAG '%1' is duplicated! (on lines %2 and %3)")
                                                    .arg(flag)
                                                    .arg(flag_to_line[flag].first)
                                                    .arg(i),
                                                CODELOC);

                // skip the FLAG line, and the FORMAT line that must come immediately after
                flag_to_line.insert(flag, QPair<qint64, qint64>(i + 2, -1));
                last_flag = flag;
                score += 1;
            }
        }
    }

    if (not last_flag.isNull())
    {
        flag_to_line[last_flag].second = nlines - flag_to_line[last_flag].first;
        last_flag = QString();
    }

    // now process all of the flag data
    score += this->processAllFlags();

    if (map["forcefield"].hasValue())
    {
        ffield = map["forcefield"].value().asA<MMDetail>();

        if (not ffield.isAmberStyle() and not ffield.isOPLS())
            throw SireError::incompatible_error(
                QObject::tr("This AmberPrm reader can only parse Amber parm files that hold molecules "
                            "that are parameterised using an Amber-style forcefield. It cannot read "
                            "molecules using the forcefield\n%1")
                    .arg(ffield.toString()),
                CODELOC);

        if (ffield.usesGeometricCombiningRules())
        {
            comb_rules = "geometric";
        }
        else
        {
            comb_rules = "arithmetic";
        }
    }
    else
    {
        // now guess the forcefield based on what we know about the potential
        ffield = MMDetail::guessFrom(comb_rules, "coulomb", "lj", 1.0 / 1.2, 0.5, "harmonic", "harmonic", "cosine");
    }

    // finally, make sure that we have been constructed sane
    this->assertSane();

    this->setScore(score);
}

/** Construct by reading from the file called 'filename' */
AmberPrm::AmberPrm(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<AmberPrm, MoleculeParser>(filename, map)
{
    this->parse(map);
}

/** Construct by reading from the contained in the passed
    set of lines */
AmberPrm::AmberPrm(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<AmberPrm, MoleculeParser>(lines, map)
{
    this->parse(map);
}

/** Internal function used to get the atom names in AtomIdx order from the
    passed AmberParams object */
QVector<QString> getAtomNames(const AmberParams &params)
{
    const int natoms = params.info().nAtoms();

    QVector<QString> atoms(natoms);

    for (int i = 0; i < natoms; ++i)
    {
        atoms[i] = params.info().name(AtomIdx(i)).value().simplified();
    }

    return atoms;
}

/** Internal function used to get the atom charges in AtomIdx order from
    the passed AmberParams object */
QVector<double> getAtomCharges(const AmberParams &params)
{
    const auto info = params.info();

    QVector<double> charges(info.nAtoms());

    const auto c = params.charges();

    for (int i = 0; i < info.nAtoms(); ++i)
    {
        charges[i] = c[info.cgAtomIdx(AtomIdx(i))].to(mod_electron) * AMBERCHARGECONV;
    }

    return charges;
}

/** Internal function used to get the atom masses in AtomIdx order from
    the passed AmberParams object */
QVector<double> getAtomMasses(const AmberParams &params)
{
    const auto info = params.info();

    QVector<double> masses(info.nAtoms());

    const auto m = params.masses();

    for (int i = 0; i < info.nAtoms(); ++i)
    {
        masses[i] = m[info.cgAtomIdx(AtomIdx(i))].to(g_per_mol);
    }

    return masses;
}

/** Internal function used to get the atom radii in AtomIdx order from
    the passed AmberParams object */
QVector<double> getAtomRadii(const AmberParams &params)
{
    const auto info = params.info();

    QVector<double> radii(info.nAtoms(), 0.0);

    const auto r = params.gbRadii();

    if (not r.isEmpty())
    {
        for (int i = 0; i < info.nAtoms(); ++i)
        {
            radii[i] = r[info.cgAtomIdx(AtomIdx(i))].to(angstrom);
        }
    }

    return radii;
}

/** Internal function used to get the atom GB screening data in AtomIdx order from
    the passed AmberParams object */
QVector<double> getAtomScreening(const AmberParams &params)
{
    const auto info = params.info();

    QVector<double> screening(info.nAtoms(), 0.0);

    const auto s = params.gbScreening();

    if (not s.isEmpty())
    {
        for (int i = 0; i < info.nAtoms(); ++i)
        {
            screening[i] = s[info.cgAtomIdx(AtomIdx(i))];
        }
    }

    return screening;
}

/** Internal function used to get the tree chain information for each atom
    in AtomIdx order from the passed AmberParams object */
QVector<QString> getAtomTreeChains(const AmberParams &params)
{
    const auto info = params.info();

    QVector<QString> treechains(info.nAtoms(), "E");

    const auto t = params.treeChains();

    if (not t.isEmpty())
    {
        for (int i = 0; i < info.nAtoms(); ++i)
        {
            treechains[i] = t[info.cgAtomIdx(AtomIdx(i))];
        }
    }

    return treechains;
}

/** Internal function used to get the amber atom types of each atom */
QVector<QString> getAmberTypes(const AmberParams &params, QHash<QString, int> &unique_names)
{
    const auto info = params.info();

    QVector<QString> ambtyps(info.nAtoms());

    const auto t = params.amberTypes();

    for (int i = 0; i < info.nAtoms(); ++i)
    {
        auto atomtype = t[info.cgAtomIdx(AtomIdx(i))];

        if (atomtype.length() > 4)
        {
            // we need to shorten the atom type name to 4 characters. Do this by
            // creating a unique name...
            int index = 0;

            if (unique_names.contains(atomtype))
            {
                index = unique_names[atomtype];
            }
            else
            {
                // find the highest index used so far
                int highest = 0;

                for (auto it = unique_names.constBegin(); it != unique_names.constEnd(); ++it)
                {
                    if (it.value() > highest)
                        highest = it.value();
                }

                highest += 1;
                unique_names.insert(atomtype, highest);

                index = highest;
            }

            if (index < 1000)
                atomtype = QString("X%1").arg(index);
            else if (index < 10000)
                atomtype = QString("%1").arg(index);
            else
                throw SireError::incompatible_error(
                    QObject::tr("Cannot write an Amber file as the number of atom types with string "
                                "lengths greater than 4 is greater than 9999"),
                    CODELOC);
        }

        ambtyps[i] = atomtype;
    }

    return ambtyps;
}

/** Internal function used to get the atom numbers */
QVector<qint64> getAtomNumbers(const AmberParams &params)
{
    const auto info = params.info();

    QVector<qint64> elements(info.nAtoms());

    const auto e = params.elements();

    for (int i = 0; i < info.nAtoms(); ++i)
    {
        elements[i] = e[info.cgAtomIdx(AtomIdx(i))].nProtons();
    }

    return elements;
}

/** Internal function used to get all of the residue data */
std::tuple<QVector<QString>, QVector<qint64>> getResidueData(const AmberParams &params, int start_idx)
{
    const auto info = params.info();

    const int nres = info.nResidues();

    if (nres == 0)
        return std::tuple<QVector<QString>, QVector<qint64>>();

    QVector<QString> res_names(nres);
    QVector<qint64> res_atoms(nres);

    for (int i = 0; i < nres; ++i)
    {
        res_names[i] = info.name(ResIdx(i)).value();

        // get the index of the first atom in this residue. As atoms are
        // contiguous (validated before here) we know that this is the
        // right pointer value for this atom
        AtomIdx atomidx = info.getAtom(ResIdx(i), 0);

        // amber uses 1-indexing
        res_atoms[i] = atomidx.value() + start_idx + 1;
    }

    return std::make_tuple(res_names, res_atoms);
}

/** Internal function used to get the excluded atoms from the
    non-bonded pairs */
QVector<QVector<qint64>> getExcludedAtoms(const AmberParams &params, int start_idx, bool use_parallel = false)
{
    const auto info = params.info();

    QVector<QVector<qint64>> excl(info.nAtoms());

    auto excl_data = excl.data();

    const auto nbpairs = params.excludedAtoms();

    if (use_parallel and info.nCutGroups() > 10)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, info.nCutGroups()), [&](const tbb::blocked_range<int> &r)
                          {
            for (int icg = r.begin(); icg < r.end(); ++icg)
            {
                const auto excluded_cg_atoms = nbpairs.excludedAtoms(CGIdx(icg));

                for (auto it = excluded_cg_atoms.constBegin();
                        it != excluded_cg_atoms.constEnd();
                        ++it)
                {
                    const int i = it.key().value();
                    const auto &excluded_atoms = it.value();

                    QVector<qint64> excl_atoms;
                    excl_atoms.reserve(excluded_atoms.count());

                    for (int j = 0; j < excluded_atoms.count(); ++j)
                    {
                        // only include i,j pairs where i < j
                        if (i < excluded_atoms[j].value())
                        {
                            // have to add 1 as Amber is 1-indexed, and add start_idx as
                            // each molecule is listed in sequence with increasing atom number
                            excl_atoms.append(excluded_atoms[j].value() + start_idx + 1);
                        }
                    }

                    if (excl_atoms.isEmpty())
                    {
                        // amber cannot have an empty excluded atoms list. In these cases
                        // it has a single atom with index 0 (0 is a null atom in amber)
                        excl_atoms.append(0);
                    }

                    excl_data[i] = excl_atoms;
                }
            } });
    }
    else if (info.nCutGroups() == 1 and info.nAtoms() <= 3)
    {
        // this is a small molecule (e.g. solvent) where all atoms
        // are excluded
        const int nats = info.nAtoms();

        if (nats == 1)
        {
            excl_data[0] = QVector<qint64>(1, 0);
        }
        else if (nats == 2)
        {
            // only add pairs where i < j
            excl_data[0] = QVector<qint64>(1, start_idx + 2);
            excl_data[1] = QVector<qint64>(1, 0);
        }
        else if (nats == 3)
        {
            // only add pairs where i < j
            excl_data[0] = QVector<qint64>(2, start_idx + 2);
            excl_data[0][1] += 1;
            excl_data[1] = QVector<qint64>(1, start_idx + 3);
            excl_data[2] = QVector<qint64>(1, 0);
        }
    }
    else
    {
        for (int icg = 0; icg < info.nCutGroups(); ++icg)
        {
            const auto excluded_cg_atoms = nbpairs.excludedAtoms(CGIdx(icg));

            for (auto it = excluded_cg_atoms.constBegin();
                 it != excluded_cg_atoms.constEnd();
                 ++it)
            {
                const int i = it.key().value();
                const auto &excluded_atoms = it.value();

                QVector<qint64> excl_atoms;
                excl_atoms.reserve(excluded_atoms.count());

                for (int j = 0; j < excluded_atoms.count(); ++j)
                {
                    // only include i,j pairs where i < j
                    if (i < excluded_atoms[j].value())
                    {
                        // have to add 1 as Amber is 1-indexed, and add start_idx as
                        // each molecule is listed in sequence with increasing atom number
                        excl_atoms.append(excluded_atoms[j].value() + start_idx + 1);
                    }
                }

                if (excl_atoms.isEmpty())
                {
                    // amber cannot have an empty excluded atoms list. In these cases
                    // it has a single atom with index 0 (0 is a null atom in amber)
                    excl_atoms.append(0);
                }

                excl_data[i] = excl_atoms;
            }
        }
    }

    return excl;
}

namespace detail
{
    /** Internal class that is used to help sort the bond arrays */
    struct Idx3
    {
    public:
        qint64 a, b, c;

        bool operator<(const Idx3 &other) const
        {
            if (a < other.a)
            {
                return true;
            }
            else if (a == other.a)
            {
                if (b < other.b)
                {
                    return true;
                }
                else if (b == other.b)
                {
                    return c < other.c;
                }
                else
                    return false;
            }
            return false;
        }
    };
} // namespace detail

/** Internal function used to get the bond information from the passed molecule */
std::tuple<QVector<qint64>, QVector<qint64>, QHash<AmberBond, qint64>> getBondData(const AmberParams &params,
                                                                                   int start_idx)
{
    QHash<AmberBond, qint64> param_to_idx;
    QVector<qint64> bonds_inc_h, bonds_exc_h;

    const auto info = params.info();
    const auto bonds = params.bonds();

    bonds_inc_h.reserve(bonds.count());
    bonds_exc_h.reserve(bonds.count());

    for (auto it = bonds.constBegin(); it != bonds.constEnd(); ++it)
    {
        // have we seen this bond parameter before?
        qint64 idx = param_to_idx.value(it.value().first, -1);

        // if not, save this to the database and increase the number of parameters
        if (idx == -1)
        {
            param_to_idx.insert(it.value().first, param_to_idx.count() + 1);
            idx = param_to_idx.count();
        }

        // is the bond marked as including hydrogen?
        if (it.value().second)
        {
            // The true atom number equals the absolute value of the number
            // divided by three, plus one (plus one as amber is 1-indexed).
            bonds_inc_h.append(3 * (info.atomIdx(it.key().atom0()).value() + start_idx));
            bonds_inc_h.append(3 * (info.atomIdx(it.key().atom1()).value() + start_idx));
            bonds_inc_h.append(idx);
        }
        else
        {
            bonds_exc_h.append(3 * (info.atomIdx(it.key().atom0()).value() + start_idx));
            bonds_exc_h.append(3 * (info.atomIdx(it.key().atom1()).value() + start_idx));
            bonds_exc_h.append(idx);
        }
    }

    // now sort the arrays so that the order is the same every time this
    // molecule is written to a file
    ::detail::Idx3 *start_it = reinterpret_cast<::detail::Idx3 *>(bonds_inc_h.data());
    ::detail::Idx3 *end_it = start_it + (bonds_inc_h.count() / 3);
    std::sort(start_it, end_it);

    start_it = reinterpret_cast<::detail::Idx3 *>(bonds_exc_h.data());
    end_it = start_it + (bonds_exc_h.count() / 3);
    std::sort(start_it, end_it);

    return std::make_tuple(bonds_inc_h, bonds_exc_h, param_to_idx);
}

namespace detail
{
    /** Internal class that is used to help sort the angle arrays */
    struct Idx4
    {
    public:
        qint64 a, b, c, d;

        bool operator<(const Idx4 &other) const
        {
            if (a < other.a)
            {
                return true;
            }
            else if (a == other.a)
            {
                if (b < other.b)
                {
                    return true;
                }
                else if (b == other.b)
                {
                    if (c < other.c)
                    {
                        return true;
                    }
                    else if (c == other.c)
                    {
                        return d < other.d;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            else
                return false;
        }
    };
} // namespace detail

/** Internal function used to get the angle information from the passed molecule */
std::tuple<QVector<qint64>, QVector<qint64>, QHash<AmberAngle, qint64>> getAngleData(const AmberParams &params,
                                                                                     int start_idx)
{
    QHash<AmberAngle, qint64> param_to_idx;
    QVector<qint64> angs_inc_h, angs_exc_h;

    const auto info = params.info();
    const auto angles = params.angles();

    angs_inc_h.reserve(angles.count());
    angs_exc_h.reserve(angles.count());

    for (auto it = angles.constBegin(); it != angles.constEnd(); ++it)
    {
        // have we seen this angle parameter before?
        qint64 idx = param_to_idx.value(it.value().first, -1);

        // if not, save this to the database and increase the number of parameters
        if (idx == -1)
        {
            param_to_idx.insert(it.value().first, param_to_idx.count() + 1);
            idx = param_to_idx.count();
        }

        // is the angle marked as including hydrogen?
        if (it.value().second)
        {
            // The true atom number equals the absolute value of the number
            // divided by three, plus one (plus one as amber is 1-indexed).
            angs_inc_h.append(3 * (info.atomIdx(it.key().atom0()).value() + start_idx));
            angs_inc_h.append(3 * (info.atomIdx(it.key().atom1()).value() + start_idx));
            angs_inc_h.append(3 * (info.atomIdx(it.key().atom2()).value() + start_idx));
            angs_inc_h.append(idx);
        }
        else
        {
            angs_exc_h.append(3 * (info.atomIdx(it.key().atom0()).value() + start_idx));
            angs_exc_h.append(3 * (info.atomIdx(it.key().atom1()).value() + start_idx));
            angs_exc_h.append(3 * (info.atomIdx(it.key().atom2()).value() + start_idx));
            angs_exc_h.append(idx);
        }
    }

    // now sort the arrays so that the order is the same every time this
    // molecule is written to a file
    ::detail::Idx4 *start_it = reinterpret_cast<::detail::Idx4 *>(angs_inc_h.data());
    ::detail::Idx4 *end_it = start_it + (angs_inc_h.count() / 4);
    std::sort(start_it, end_it);

    start_it = reinterpret_cast<::detail::Idx4 *>(angs_exc_h.data());
    end_it = start_it + (angs_exc_h.count() / 4);
    std::sort(start_it, end_it);

    return std::make_tuple(angs_inc_h, angs_exc_h, param_to_idx);
}

namespace detail
{
    /** Internal class that is used to help sort the dihedral arrays */
    struct Idx5
    {
    public:
        qint64 a, b, c, d, e;

        bool operator<(const Idx5 &other) const
        {
            qint64 c_cmp = std::abs(c);
            qint64 other_c_cmp = std::abs(other.c);
            qint64 d_cmp = std::abs(d);
            qint64 other_d_cmp = std::abs(other.d);
            if (a < other.a)
            {
                return true;
            }
            else if (a == other.a)
            {
                if (b < other.b)
                {
                    return true;
                }
                else if (b == other.b)
                {
                    if (c_cmp < other_c_cmp)
                    {
                        return true;
                    }
                    else if (c_cmp == other_c_cmp)
                    {
                        if (d_cmp < other_d_cmp)
                        {
                            return true;
                        }
                        else if (d_cmp == other_d_cmp)
                        {
                            return e < other.e;
                        }
                    }
                }
            }
            return false;
        }
    };
} // namespace detail

/** Internal function used to get the dihedral information from the passed molecule */
std::tuple<QVector<qint64>, QVector<qint64>, QHash<AmberNBDihPart, qint64>> getDihedralData(const AmberParams &params,
                                                                                            int start_idx)
{
    QHash<AmberNBDihPart, qint64> param_to_idx;
    QHash<qint64, AmberNBDihPart> idx_to_param;
    QVector<qint64> dihs_inc_h, dihs_exc_h;

    const auto info = params.info();
    const auto dihedrals = params.dihedrals();
    const auto impropers = params.impropers();

    dihs_inc_h.reserve(10 * (dihedrals.count() + impropers.count()));
    dihs_exc_h.reserve(10 * (dihedrals.count() + impropers.count()));

    for (auto it = dihedrals.constBegin(); it != dihedrals.constEnd(); ++it)
    {
        // get the NB14 scaling factor for this dihedral. This should be set
        // for only the first dihedral found that has this pair, otherwise
        // we would risk Amber double-counting this parameter
        const BondID nb14pair = params.convert(BondID(it.key().atom0(), it.key().atom3()));
        const AmberNB14 nb14 = params.nb14s().value(nb14pair);

        // extract the individual dihedral terms from the dihedral,
        // combine them with the NB14 scaling factors, and then
        // create a database of them and get their ID
        QList<qint64> idxs;
        QList<qint64> ignored;

        if (it.value().first.terms().isEmpty())
        {
            // all dihedrals must have some terms, or we risk losing
            //  NB14 scale factors
            throw SireError::program_bug(QObject::tr("How can a dihedral have no terms? %1").arg(it.key().toString()),
                                         CODELOC);
        }

        for (const auto &term : it.value().first.terms())
        {
            AmberNBDihPart nbterm(term, nb14);

            if (nb14.cscl() == 0 and nb14.ljscl() == 0)
            {
                ignored.append(-1);
            }
            else
            {
                // this is a non-zero coulomb and LJ scale,
                // so this interaction should not be ignored
                ignored.append(1);
            }

            qint64 idx = param_to_idx.value(nbterm, -1);

            if (idx == -1)
            {
                idx = param_to_idx.count() + 1;
                param_to_idx.insert(nbterm, idx);
                idx_to_param.insert(idx, nbterm);
            }

            idxs.append(idx);
        }

        // The true atom number equals the absolute value of the number
        // divided by three, plus one (plus one as amber is 1-indexed).
        qint64 atom0 = 3 * (info.atomIdx(it.key().atom0()).value() + start_idx);
        qint64 atom1 = 3 * (info.atomIdx(it.key().atom1()).value() + start_idx);
        qint64 atom2 = 3 * (info.atomIdx(it.key().atom2()).value() + start_idx);
        qint64 atom3 = 3 * (info.atomIdx(it.key().atom3()).value() + start_idx);

        if (atom2 == 0 or atom3 == 0)
        {
            // we can't record a negative zero, so we need to invert the
            // indicies
            std::swap(atom0, atom3);
            std::swap(atom1, atom2);
        }

        // is the dihedral marked as including hydrogen?
        if (it.value().second)
        {
            for (int i = 0; i < idxs.count(); ++i)
            {
                dihs_inc_h.append(atom0);
                dihs_inc_h.append(atom1);
                dihs_inc_h.append(ignored[i] * atom2);
                dihs_inc_h.append(atom3);
                dihs_inc_h.append(idxs[i]);
            }
        }
        else
        {
            for (int i = 0; i < idxs.count(); ++i)
            {
                dihs_exc_h.append(atom0);
                dihs_exc_h.append(atom1);
                dihs_exc_h.append(ignored[i] * atom2);
                dihs_exc_h.append(atom3);
                dihs_exc_h.append(idxs[i]);
            }
        }
    }

    for (auto it = impropers.constBegin(); it != impropers.constEnd(); ++it)
    {
        // extract the individual dihedral terms from the dihedral parameters,
        // combine them with the 0,0 scaling factor, and then
        // create a database of them and get their ID
        QList<qint64> idxs;

        for (const auto &term : it.value().first.terms())
        {
            AmberNBDihPart nbterm(term, AmberNB14(0, 0));

            qint64 idx = param_to_idx.value(nbterm, -1);

            if (idx == -1)
            {
                idx = param_to_idx.count() + 1;
                param_to_idx.insert(nbterm, idx);
                idx_to_param.insert(idx, nbterm);
            }

            idxs.append(idx);
        }

        // The true atom number equals the absolute value of the number
        // divided by three, plus one (plus one as amber is 1-indexed).
        qint64 atom0 = 3 * (info.atomIdx(it.key().atom0()).value() + start_idx);
        qint64 atom1 = 3 * (info.atomIdx(it.key().atom1()).value() + start_idx);
        qint64 atom2 = 3 * (info.atomIdx(it.key().atom2()).value() + start_idx);
        qint64 atom3 = 3 * (info.atomIdx(it.key().atom3()).value() + start_idx);

        // is the improper marked as including hydrogen?
        if (it.value().second)
        {
            for (const auto idx : idxs)
            {
                dihs_inc_h.append(atom0);
                dihs_inc_h.append(atom1);
                dihs_inc_h.append(-atom2); // 1-4 interactions of impropers are ignored
                dihs_inc_h.append(-atom3); // negative to indicate this is an improper
                dihs_inc_h.append(idx);
            }
        }
        else
        {
            for (const auto idx : idxs)
            {
                dihs_exc_h.append(atom0);
                dihs_exc_h.append(atom1);
                dihs_exc_h.append(-atom2); // 1-4 interactions of impropers are ignored
                dihs_exc_h.append(-atom3); // negative to indicate this is an improper
                dihs_exc_h.append(idx);
            }
        }
    }

    // now sort the arrays so that the order is the same every time this
    // molecule is written to a file
    ::detail::Idx5 *start_it = reinterpret_cast<::detail::Idx5 *>(dihs_inc_h.data());
    ::detail::Idx5 *end_it = start_it + (dihs_inc_h.count() / 5);
    std::sort(start_it, end_it);

    start_it = reinterpret_cast<::detail::Idx5 *>(dihs_exc_h.data());
    end_it = start_it + (dihs_exc_h.count() / 5);
    std::sort(start_it, end_it);

    return std::make_tuple(dihs_inc_h, dihs_exc_h, param_to_idx);
}

namespace detail
{
    /** Internal class that is used to help sort the cmap arrays */
    struct Idx6
    {
    public:
        qint64 a, b, c, d, e, f;

        bool operator<(const Idx6 &other) const
        {
            qint64 c_cmp = std::abs(c);
            qint64 other_c_cmp = std::abs(other.c);
            qint64 d_cmp = std::abs(d);
            qint64 other_d_cmp = std::abs(other.d);
            qint64 e_cmp = std::abs(e);
            qint64 other_e_cmp = std::abs(other.e);
            if (a < other.a)
            {
                return true;
            }
            else if (a == other.a)
            {
                if (b < other.b)
                {
                    return true;
                }
                else if (b == other.b)
                {
                    if (c_cmp < other_c_cmp)
                    {
                        return true;
                    }
                    else if (c_cmp == other_c_cmp)
                    {
                        if (d_cmp < other_d_cmp)
                        {
                            return true;
                        }
                        else if (d_cmp == other_d_cmp)
                        {
                            if (e_cmp < other_e_cmp)
                            {
                                return true;
                            }
                            else if (e_cmp == other_e_cmp)
                            {
                                return f < other.f;
                            }
                        }
                    }
                }
            }
            return false;
        }
    };
} // namespace detail

/** Internal function used to get the cmap information from the passed molecule */
std::tuple<QVector<qint64>, QHash<CMAPParameter, qint64>> getCMAPData(const AmberParams &params,
                                                                      int start_idx)
{
    QHash<CMAPParameter, qint64> param_to_idx;
    QVector<qint64> cmap_idxs;

    start_idx /= 3;

    const auto info = params.info();
    const auto cmaps = params.cmapFunctions().parameters();

    if (cmaps.isEmpty())
    {
        return std::make_tuple(cmap_idxs, param_to_idx);
    }

    cmap_idxs.reserve(cmaps.count());

    for (auto it = cmaps.constBegin(); it != cmaps.constEnd(); ++it)
    {
        // have we seen this cmap parameter before?
        qint64 idx = param_to_idx.value(it->parameter(), -1);

        // if not, save this to the database and increase the number of parameters
        if (idx == -1)
        {
            param_to_idx.insert(it->parameter(), param_to_idx.count() + 1);
            idx = param_to_idx.count();
        }

        // need to add 1 as Amber is 1-indexed
        cmap_idxs.append(info.atomIdx(it->atom0()).value() + 1 + start_idx);
        cmap_idxs.append(info.atomIdx(it->atom1()).value() + 1 + start_idx);
        cmap_idxs.append(info.atomIdx(it->atom2()).value() + 1 + start_idx);
        cmap_idxs.append(info.atomIdx(it->atom3()).value() + 1 + start_idx);
        cmap_idxs.append(info.atomIdx(it->atom4()).value() + 1 + start_idx);
        cmap_idxs.append(idx);
    }

    // now sort the arrays so that the order is the same every time this
    // molecule is written to a file
    ::detail::Idx6 *start_it = reinterpret_cast<::detail::Idx6 *>(cmap_idxs.data());
    ::detail::Idx6 *end_it = start_it + (cmap_idxs.count() / 6);
    std::sort(start_it, end_it);

    return std::make_tuple(cmap_idxs, param_to_idx);
}

/** Internal function that converts the passed list of parameters into a list
    of text lines */
QStringList toLines(const QVector<AmberParams> &params, const Space &space, int num_dummies,
                    bool use_geometric_rules, bool use_parallel = true,
                    QStringList *all_errors = 0)
{
    if (params.isEmpty())
        return QStringList();

    tbb::spin_mutex error_mutex;

    // before we start, validate that all of the parameters are ok
    bool ok = true;

    const auto params_data = params.constData();

    if (use_parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> &r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                QStringList errors = params_data[i].validate();

                if (not errors.isEmpty())
                {
                    tbb::spin_mutex::scoped_lock locker(error_mutex);
                    ok = false;

                    if (all_errors)
                    {
                        all_errors->append(QObject::tr(" == Errors validating parameters for molecule %1 of %2 == ")
                                               .arg(i + 1)
                                               .arg(params.count()));
                        (*all_errors) += errors;
                    }
                }
            } });
    }
    else
    {
        for (int i = 0; i < params.count(); ++i)
        {
            QStringList errors = params_data[i].validate();

            if (not errors.isEmpty())
            {
                ok = false;

                if (all_errors)
                {
                    all_errors->append(QObject::tr(" == Errors validating parameters for molecule %1 of %2 == ")
                                           .arg(i + 1)
                                           .arg(params.count()));
                    (*all_errors) += errors;
                }
            }
        }
    }

    if (not ok)
    {
        return QStringList();
    }

    // now set up the 'pointers' array, which forms the header of the output
    // file and gives information about the number of atoms etc.
    QVector<qint64> pointers(33, 0);
    auto pointers_data = pointers.data();

    // first, count up the number of atoms, and get the start index of each
    // atom (the atoms count sequentially upwards)
    QVector<qint64> _tmp_natoms_per_molecule(params.count(), 0);
    QVector<qint64> _tmp_atom_start_index(params.count(), 0);
    {
        int natoms = 0;

        _tmp_natoms_per_molecule[0] = params_data[0].info().nAtoms();
        natoms = _tmp_natoms_per_molecule[0];

        for (int i = 1; i < params.count(); ++i)
        {
            const int nats = params_data[i].info().nAtoms();
            _tmp_natoms_per_molecule[i] = nats;
            _tmp_atom_start_index[i] = natoms;

            natoms += nats;
        }

        pointers[0] = natoms;
    }

    // save this data as const vectors, so that we always use the thread-safe
    // const functions to access the data
    const QVector<qint64> natoms_per_molecule = _tmp_natoms_per_molecule;
    const QVector<qint64> atom_start_index = _tmp_atom_start_index;

    const auto atom_start_index_data = atom_start_index.constData();

    // save the number of atoms to be written to the file
    const int total_natoms = pointers[0];

    // now see if there is a periodic box
    const bool has_periodic_box = space.isPeriodic();

    if (has_periodic_box)
    {
        // Orthorhombic box.
        if (space.isA<PeriodicBox>())
        {
            pointers[27] = 1;
        }
        // General triclinic box.
        else if (space.isA<TriclinicBox>())
        {
            pointers[27] = 3;
        }
    }

    // here is the number of solvent molecules, and the index of the last
    // solute residue. These will be updated by the 'getAllAtomNames' function
    int nsolvents = 0;
    int last_solute_residue = 0;

    // function used to generate the text for the atom names
    auto getAllAtomNames = [&]()
    {
        QVector<QVector<QString>> atom_names(params.count());
        auto atom_names_data = atom_names.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    atom_names_data[i] = getAtomNames(params_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                atom_names_data[i] = getAtomNames(params_data[i]);
            }
        }

        // now loop from the last molecule backwards to find the first solvent
        // molecule. Solvent molecules are those that have the same number
        // of atoms as the last molecule, and that have the same atom names
        {
            const auto last_molinfo = params_data[params.count() - 1].info();

            int first_solvent = params.count() - 1;

            for (int i = params.count() - 2; i >= 0; --i)
            {
                const auto molinfo = params_data[i].info();

                // compare the number of atoms in this molecule to that of
                // the last molecule - if different, then this is not the solvent
                if (last_molinfo.nAtoms() != molinfo.nAtoms())
                {
                    break;
                }

                // compare the atom names - these must be the same in the same order
                for (AtomIdx j(0); j < molinfo.nAtoms(); ++j)
                {
                    if (last_molinfo.name(j) != molinfo.name(j))
                    {
                        break;
                    }
                }

                first_solvent = i;
            }

            nsolvents = params.count() - first_solvent;

            // now count up the number of residues in the "solute" molecules
            for (int i = 0; i < first_solvent; ++i)
            {
                last_solute_residue += params_data[i].info().nResidues();
            }
        }

        if (all_errors)
        {
            QStringList errors;

            QStringList lines = writeStringData(collapse(atom_names), AmberFormat(AmberPrm::STRING, 20, 4), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing atom name data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeStringData(collapse(atom_names), AmberFormat(AmberPrm::STRING, 20, 4));
    };

    // function used to generate the text for all atom charges
    auto getAllAtomCharges = [&]()
    {
        QVector<QVector<double>> atom_charges(params.count());
        auto atom_charges_data = atom_charges.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    atom_charges_data[i] = getAtomCharges(params_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                atom_charges_data[i] = getAtomCharges(params_data[i]);
            }
        }

        if (all_errors)
        {
            QStringList errors;

            QStringList lines = writeFloatData(collapse(atom_charges), AmberFormat(AmberPrm::FLOAT, 5, 16, 8), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing atom charge data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeFloatData(collapse(atom_charges), AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
    };

    // function used to generate the text for all atom numbers
    auto getAllAtomNumbers = [&]()
    {
        QVector<QVector<qint64>> atom_numbers(params.count());
        auto atom_numbers_data = atom_numbers.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    atom_numbers_data[i] = getAtomNumbers(params_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                atom_numbers_data[i] = getAtomNumbers(params_data[i]);
            }
        }

        if (all_errors)
        {
            QStringList errors;

            QStringList lines = writeIntData(collapse(atom_numbers), AmberFormat(AmberPrm::INTEGER, 10, 8), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing atom number data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeIntData(collapse(atom_numbers), AmberFormat(AmberPrm::INTEGER, 10, 8));
    };

    // function used to generate the text for all atom masses
    auto getAllAtomMasses = [&]()
    {
        QVector<QVector<double>> atom_masses(params.count());
        auto atom_masses_data = atom_masses.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    atom_masses_data[i] = getAtomMasses(params_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                atom_masses_data[i] = getAtomMasses(params_data[i]);
            }
        }

        if (all_errors)
        {
            QStringList errors;

            QStringList lines = writeFloatData(collapse(atom_masses), AmberFormat(AmberPrm::FLOAT, 5, 16, 8), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing atom mass data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeFloatData(collapse(atom_masses), AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
    };

    // function used to generate the text for all atom Born radii
    auto getAllRadii = [&]()
    {
        QVector<QVector<double>> atom_radii(params.count());
        auto atom_radii_data = atom_radii.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    atom_radii_data[i] = getAtomRadii(params_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                atom_radii_data[i] = getAtomRadii(params_data[i]);
            }
        }

        if (all_errors)
        {
            QStringList errors;

            QStringList lines = writeFloatData(collapse(atom_radii), AmberFormat(AmberPrm::FLOAT, 5, 16, 8), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing atom radii data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeFloatData(collapse(atom_radii), AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
    };

    // function used to generate the text for all atom Generalised Born screening parameters
    auto getAllScreening = [&]()
    {
        QVector<QVector<double>> atom_screening(params.count());
        auto atom_screening_data = atom_screening.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    atom_screening_data[i] = getAtomScreening(params_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                atom_screening_data[i] = getAtomScreening(params_data[i]);
            }
        }

        if (all_errors)
        {
            QStringList errors;

            QStringList lines =
                writeFloatData(collapse(atom_screening), AmberFormat(AmberPrm::FLOAT, 5, 16, 8), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing atom screening data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeFloatData(collapse(atom_screening), AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
    };

    // function used to generate the text for all atom tree chain information
    auto getAllTreeChains = [&]()
    {
        QVector<QVector<QString>> atom_treechains(params.count());
        auto atom_treechains_data = atom_treechains.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    atom_treechains_data[i] = getAtomTreeChains(params_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                atom_treechains_data[i] = getAtomTreeChains(params_data[i]);
            }
        }

        if (all_errors)
        {
            QStringList errors;

            QStringList lines =
                writeStringData(collapse(atom_treechains), AmberFormat(AmberPrm::STRING, 20, 4), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing atom treechain data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeStringData(collapse(atom_treechains), AmberFormat(AmberPrm::STRING, 20, 4));
    };

    // function used to get all of the Amber atom types
    auto getAllAmberTypes = [&]()
    {
        QVector<QVector<QString>> amber_types(params.count());
        auto amber_types_data = amber_types.data();

        QHash<QString, int> unique_types;  // Unique types in a molecule.
        QSet<QString> system_unique_types; // Unique types in the system.

        for (int i = 0; i < params.count(); ++i)
        {
            amber_types_data[i] = getAmberTypes(params_data[i], unique_types);

            // Store the unique types for the entire system.
            for (const auto &type : amber_types_data[i])
            {
                // We haven't seen this type in any molecule.
                if (not system_unique_types.contains(type))
                {
                    system_unique_types.insert(type);
                }
            }
        }

        // Record the NATYP flag.
        const int ntypes_system = system_unique_types.count();
        pointers_data[18] = ntypes_system;

        if (all_errors)
        {
            QStringList errors;

            QStringList lines = writeStringData(collapse(amber_types), AmberFormat(AmberPrm::STRING, 20, 4), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing amber type data! =="));
                *all_errors += errors;
            }

            return lines;
        }
        else
            return writeStringData(collapse(amber_types), AmberFormat(AmberPrm::STRING, 20, 4));
    };

    // function used to get all of the LJ parameters and atom types
    auto getAllAtomTypes = [&]()
    {
        // we need to go through each atom and work out if the atom type
        // is unique
        QVector<std::pair<LJParameter, QList<LJException>>> ljparams;
        QVector<QVector<qint64>> atom_types(params.count());
        auto atom_types_data = atom_types.data();

        bool has_ccoeff = false;
        bool has_ljexceptions = false;

        for (int i = 0; i < params.count(); ++i)
        {
            const auto info = params_data[i].info();
            const auto ljs = params_data[i].ljs();

            const bool mol_has_lj_exceptions = ljs.hasExceptions();

            QVector<qint64> mol_atom_types(info.nAtoms());
            auto mol_atom_types_data = mol_atom_types.data();

            for (int j = 0; j < info.nAtoms(); ++j)
            {
                const LJParameter lj = ljs[info.cgAtomIdx(AtomIdx(j))];

                QList<LJException> ljexceptions;

                if (mol_has_lj_exceptions)
                {
                    has_ljexceptions = true;
                    ljexceptions = ljs.getExceptions(j);

                    if (not has_ccoeff)
                    {
                        for (const auto &ljexception : ljexceptions)
                        {
                            if (ljexception.value().hasC())
                            {
                                has_ccoeff = true;
                                break;
                            }
                        }
                    }
                }

                auto lj_with_exceptions = std::make_pair(lj, ljexceptions);

                int idx = ljparams.indexOf(lj_with_exceptions);

                if (idx == -1)
                {
                    ljparams.append(lj_with_exceptions);
                    mol_atom_types_data[j] = ljparams.count();
                }
                else
                {
                    mol_atom_types_data[j] = idx + 1;
                }
            }

            atom_types_data[i] = mol_atom_types;
        }

        // We now have all of the atom types - create the acoeff and bcoeff arrays.
        const int ntypes = ljparams.count();

        pointers_data[1] = ntypes;

        QVector<qint64> ico(ntypes * ntypes);
        QVector<double> cn1((ntypes * (ntypes + 1)) / 2);
        QVector<double> cn2((ntypes * (ntypes + 1)) / 2);
        QVector<double> cn3;

        if (has_ccoeff)
        {
            cn3 = QVector<double>((ntypes * (ntypes + 1)) / 2, 0.0);
        }

        auto ico_data = ico.data();
        auto cn1_data = cn1.data();
        auto cn2_data = cn2.data();

        const auto ljparams_data = ljparams.constData();

        int lj_idx = 0;

        for (int i = 0; i < ntypes; ++i)
        {
            const auto &lj0 = std::get<0>(ljparams_data[i]);

            // now we need to (wastefully) save the A/B parameters for
            // all mixed LJ parameters i+j
            for (int j = 0; j < i; ++j)
            {
                LJParameter lj01 = lj0.combineArithmetic(std::get<0>(ljparams.constData()[j]));

                if (use_geometric_rules)
                    lj01 = lj0.combineGeometric(std::get<0>(ljparams.constData()[j]));

                cn1_data[lj_idx] = lj01.A();
                cn2_data[lj_idx] = lj01.B();

                lj_idx += 1;

                // now save the index for i,j and j,i into the ico array
                ico_data[i * ntypes + j] = lj_idx;
                ico_data[j * ntypes + i] = lj_idx;
            }

            // now save the A/B parameters for the ith parameter type
            cn1_data[lj_idx] = lj0.A();
            cn2_data[lj_idx] = lj0.B();

            lj_idx += 1;
            ico_data[i * ntypes + i] = lj_idx;
        }

        if (has_ljexceptions)
        {
            int check_lj_idx = lj_idx;

            // go through every pair of parameters and see if they have
            // an exception - if they do, then update the arrays to have
            // the exception values
            int lj_idx = 0;

            for (int i = 0; i < ntypes; ++i)
            {
                const auto &lj_exceptions0 = std::get<1>(ljparams_data[i]);

                if (lj_exceptions0.isEmpty())
                {
                    lj_idx += (i + 1);
                    continue;
                }

                for (int j = 0; j < i; ++j)
                {
                    const auto &lj_exceptions1 = std::get<1>(ljparams_data[j]);

                    if (lj_exceptions1.isEmpty())
                    {
                        lj_idx += 1;
                        continue;
                    }

                    // see if there is a match - use the first match
                    bool found_match = false;

                    for (const auto &lj_exception0 : lj_exceptions0)
                    {
                        for (const auto &lj_exception1 : lj_exceptions1)
                        {
                            if (lj_exception0.pairsWith(lj_exception1))
                            {
                                // we have a match - update the arrays
                                found_match = true;
                                const auto &lj = lj_exception0.value();

                                cn1_data[lj_idx] = lj.A();
                                cn2_data[lj_idx] = lj.B();

                                if (lj_exception0.value().hasC())
                                {
                                    cn3[lj_idx] = lj.C();
                                }

                                break;
                            }
                        }

                        if (found_match)
                            break;
                    }

                    lj_idx += 1;
                }

                // now check for the self-exception
                bool found_match = false;
                for (const auto &lj_exception0 : lj_exceptions0)
                {
                    for (const auto &lj_exception1 : lj_exceptions0)
                    {
                        if (lj_exception0.pairsWith(lj_exception1))
                        {
                            // we have a match - update the arrays
                            found_match = true;
                            const auto &lj = lj_exception0.value();

                            cn1_data[lj_idx] = lj.A();
                            cn2_data[lj_idx] = lj.B();

                            if (lj_exception0.value().hasC())
                            {
                                cn3[lj_idx] = lj.C();
                            }

                            break;
                        }
                    }

                    if (found_match)
                        break;
                }

                lj_idx += 1;
            }

            if (lj_idx != check_lj_idx)
            {
                throw SireError::program_bug(QObject::tr("The number of LJ parameters with exceptions (%1) does not match the number of LJ parameters (%2)")
                                                 .arg(lj_idx)
                                                 .arg(check_lj_idx),
                                             CODELOC);
            }
        }

        QStringList cn3_strings;

        if (has_ccoeff)
        {
            cn3_strings = writeFloatData(cn3, AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
        }

        // now return all of the arrays
        return std::make_tuple(writeIntData(collapse(atom_types), AmberFormat(AmberPrm::INTEGER, 10, 8)),
                               writeIntData(ico, AmberFormat(AmberPrm::INTEGER, 10, 8)),
                               writeFloatData(cn1, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeFloatData(cn2, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               cn3_strings);
    };

    // function used to generate the text for the excluded atoms
    auto getAllExcludedAtoms = [&]()
    {
        const int nmols = params.count();

        QVector<QVector<qint64>> all_excluded_atoms(nmols);
        QVector<QVector<qint64>> num_excluded_atoms(nmols);

        auto all_excluded_atoms_data = all_excluded_atoms.data();
        auto num_excluded_atoms_data = num_excluded_atoms.data();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    const auto excluded_atoms = getExcludedAtoms(params_data[i], atom_start_index_data[i], true);

                    QVector<qint64> mol_all_excluded_atoms;
                    QVector<qint64> mol_num_excluded_atoms;

                    for (const auto &e : excluded_atoms)
                    {
                        mol_num_excluded_atoms.append(e.count());
                        mol_all_excluded_atoms += e;
                    }

                    all_excluded_atoms_data[i] = mol_all_excluded_atoms;
                    num_excluded_atoms_data[i] = mol_num_excluded_atoms;
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                const auto excluded_atoms = getExcludedAtoms(params_data[i], atom_start_index_data[i], false);

                QVector<qint64> mol_all_excluded_atoms;
                QVector<qint64> mol_num_excluded_atoms;

                for (const auto &e : excluded_atoms)
                {
                    mol_num_excluded_atoms.append(e.count());
                    mol_all_excluded_atoms += e;
                }

                all_excluded_atoms_data[i] = mol_all_excluded_atoms;
                num_excluded_atoms_data[i] = mol_num_excluded_atoms;
            }
        }

        QStringList excluded_atoms_lines, num_excluded_lines;

        if (all_errors)
        {
            QStringList errors;

            excluded_atoms_lines = writeIntData(collapse(all_excluded_atoms), AmberFormat(AmberPrm::INTEGER, 10, 8), &errors);

            num_excluded_lines = writeIntData(collapse(num_excluded_atoms), AmberFormat(AmberPrm::INTEGER, 10, 8), &errors);

            if (not errors.isEmpty())
            {
                tbb::spin_mutex::scoped_lock locker(error_mutex);
                all_errors->append(QObject::tr("== Errors writing excluded atom data! =="));
                *all_errors += errors;
            }
        }
        else
        {
            excluded_atoms_lines = writeIntData(collapse(all_excluded_atoms), AmberFormat(AmberPrm::INTEGER, 10, 8));
            num_excluded_lines = writeIntData(collapse(num_excluded_atoms), AmberFormat(AmberPrm::INTEGER, 10, 8));
        }

        int nexcl = 0;

        for (const auto &e : all_excluded_atoms)
        {
            nexcl += e.count();
        }

        return std::make_tuple(num_excluded_lines, excluded_atoms_lines, nexcl);
    };

    // function used to get all of the info for all of the residues
    auto getAllResidueInfo = [&]()
    {
        QVector<std::tuple<QVector<QString>, QVector<qint64>>> all_res_data(params.count());
        auto all_res_data_data = all_res_data.data();

        // get all of the residue info - this is the name and the atom range
        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    all_res_data_data[i] = getResidueData(params_data[i], atom_start_index_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                all_res_data_data[i] = getResidueData(params_data[i], atom_start_index_data[i]);
            }
        }

        auto res_data = collapseTuples(all_res_data);

        // get the number of residues
        int nres = std::get<0>(res_data).count();

        // now find out the maximum number of atoms in a residue
        int maxres = 0;

        for (int i = 1; i < nres; ++i)
        {
            int nats = std::get<1>(res_data)[i] - std::get<1>(res_data)[i - 1];

            if (nats > maxres)
                maxres = nats;
        }

        int nats = total_natoms - std::get<1>(res_data)[nres - 1] + 1;

        if (nats > maxres)
            maxres = nats;

        return std::make_tuple(writeStringData(std::get<0>(res_data), AmberFormat(AmberPrm::STRING, 20, 4)),
                               writeIntData(std::get<1>(res_data), AmberFormat(AmberPrm::INTEGER, 10, 8)), nres,
                               maxres);
    };

    // function used to get all of the bond information
    auto getAllBonds = [&]()
    {
        QVector<std::tuple<QVector<qint64>, QVector<qint64>, QHash<AmberBond, qint64>>> all_bond_data(params.count());
        auto all_bond_data_data = all_bond_data.data();

        // get all of the bond information from each molecule
        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    all_bond_data_data[i] = getBondData(params_data[i], atom_start_index_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                all_bond_data_data[i] = getBondData(params_data[i], atom_start_index_data[i]);
            }
        }

        // count up the number of bond parameters to help reserve memory space
        int nbonds_inc_h = 0;
        int nbonds_exc_h = 0;
        int nparams = 0;
        for (int i = 0; i < params.count(); ++i)
        {
            nbonds_inc_h += std::get<0>(all_bond_data[i]).count();
            nbonds_exc_h += std::get<1>(all_bond_data[i]).count();
            nparams += std::get<2>(all_bond_data[i]).count();
        }

        QVector<qint64> all_bonds_inc_h, all_bonds_exc_h;
        all_bonds_inc_h.reserve(3 * nbonds_inc_h);
        all_bonds_exc_h.reserve(3 * nbonds_exc_h);

        // we now need to go through the bonds and remove duplicated bond parameters
        QVector<double> k_data, r0_data;
        QHash<AmberBond, int> all_bond_to_idx;
        all_bond_to_idx.reserve(nparams);

        // combine all parameters into a single set, with a unique index
        for (int i = 0; i < params.count(); ++i)
        {
            const auto bond_to_idx = std::get<2>(all_bond_data_data[i]);

            // first, find all unique bonds
            for (auto it = bond_to_idx.constBegin(); it != bond_to_idx.constEnd(); ++it)
            {
                if (not all_bond_to_idx.contains(it.key()))
                {
                    all_bond_to_idx.insert(it.key(), -1);
                }
            }
        }

        // now, sort the list so that there is a consistent order of
        // parameters every time this output file is written
        {
            auto all_bonds = all_bond_to_idx.keys();
            std::sort(all_bonds.begin(), all_bonds.end());

            k_data = QVector<double>(all_bonds.count());
            r0_data = QVector<double>(all_bonds.count());

            int i = 0;
            for (auto bond : all_bonds)
            {
                k_data[i] = bond.k();
                r0_data[i] = bond.r0();

                all_bond_to_idx[bond] = i + 1;
                i += 1;
            }
        }

        for (int i = 0; i < params.count(); ++i)
        {
            // create a hash to map local parameter indicies to global parameter indicies
            QHash<qint64, qint64> idx_to_idx;

            // get all of the bond parameters
            const auto bond_to_idx = std::get<2>(all_bond_data_data[i]);

            // for each one, get the global parameter index, save the parameters to
            // the k and r0 arrays, and update the mapping from local to global indicies
            for (auto it = bond_to_idx.constBegin(); it != bond_to_idx.constEnd(); ++it)
            {
                idx_to_idx.insert(it.value(), all_bond_to_idx.value(it.key(), -1));
            }

            // now run through the bonds updating their local indicies to
            // global indicies
            const auto bonds_inc_h = std::get<0>(all_bond_data_data[i]);
            const auto bonds_exc_h = std::get<1>(all_bond_data_data[i]);

            for (int j = 0; j < bonds_inc_h.count(); j += 3)
            {
                all_bonds_inc_h.append(bonds_inc_h[j]);
                all_bonds_inc_h.append(bonds_inc_h[j + 1]);
                all_bonds_inc_h.append(idx_to_idx.value(bonds_inc_h[j + 2]));
            }

            for (int j = 0; j < bonds_exc_h.count(); j += 3)
            {
                all_bonds_exc_h.append(bonds_exc_h[j]);
                all_bonds_exc_h.append(bonds_exc_h[j + 1]);
                all_bonds_exc_h.append(idx_to_idx.value(bonds_exc_h[j + 2]));
            }
        }

        int nbonh = all_bonds_inc_h.count() / 3; // number of bonds containing hydrogen
        int mbona = all_bonds_exc_h.count() / 3; // number of bonds not containing hydrogen
        int nbona = mbona;                       // mbona + number of constraint bonds
        int numbnd = k_data.count();             // number of unique bond types

        return std::make_tuple(writeFloatData(k_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeFloatData(r0_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeIntData(all_bonds_inc_h, AmberFormat(AmberPrm::INTEGER, 10, 8)),
                               writeIntData(all_bonds_exc_h, AmberFormat(AmberPrm::INTEGER, 10, 8)), nbonh, mbona,
                               nbona, numbnd);
    };

    // function used to get all of the angle information
    auto getAllAngles = [&]()
    {
        QVector<std::tuple<QVector<qint64>, QVector<qint64>, QHash<AmberAngle, qint64>>> all_ang_data(params.count());
        auto all_ang_data_data = all_ang_data.data();

        // get all of the angle information from each molecule
        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    all_ang_data_data[i] = getAngleData(params_data[i], atom_start_index_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                all_ang_data_data[i] = getAngleData(params_data[i], atom_start_index_data[i]);
            }
        }

        // count up the number of angle parameters to help reserve memory space
        int nangs_inc_h = 0;
        int nangs_exc_h = 0;
        int nparams = 0;
        for (int i = 0; i < params.count(); ++i)
        {
            nangs_inc_h += std::get<0>(all_ang_data[i]).count();
            nangs_exc_h += std::get<1>(all_ang_data[i]).count();
            nparams += std::get<2>(all_ang_data[i]).count();
        }

        QVector<qint64> all_angs_inc_h, all_angs_exc_h;
        all_angs_inc_h.reserve(4 * nangs_inc_h);
        all_angs_exc_h.reserve(4 * nangs_exc_h);

        // we now need to go through the angles and remove duplicated angle parameters
        QVector<double> k_data, t0_data;
        QHash<AmberAngle, int> all_ang_to_idx;
        all_ang_to_idx.reserve(nparams);

        // combine all parameters into a single set, with a unique index
        for (int i = 0; i < params.count(); ++i)
        {
            const auto ang_to_idx = std::get<2>(all_ang_data_data[i]);

            // first, find all unique angles
            for (auto it = ang_to_idx.constBegin(); it != ang_to_idx.constEnd(); ++it)
            {
                if (not all_ang_to_idx.contains(it.key()))
                {
                    all_ang_to_idx.insert(it.key(), -1);
                }
            }
        }

        // now, sort the list so that there is a consistent order of
        // parameters every time this output file is written
        {
            auto all_angs = all_ang_to_idx.keys();
            std::sort(all_angs.begin(), all_angs.end());

            k_data = QVector<double>(all_angs.count());
            t0_data = QVector<double>(all_angs.count());

            int i = 0;
            for (auto ang : all_angs)
            {
                k_data[i] = ang.k();
                t0_data[i] = ang.theta0();

                all_ang_to_idx[ang] = i + 1;
                i += 1;
            }
        }

        for (int i = 0; i < params.count(); ++i)
        {
            // create a hash to map local parameter indicies to global parameter indicies
            QHash<qint64, qint64> idx_to_idx;

            // get all of the angle parameters
            const auto ang_to_idx = std::get<2>(all_ang_data_data[i]);

            // for each one, get the global parameter index, save the parameters to
            // the k and t0 arrays, and update the mapping from local to global indicies
            for (auto it = ang_to_idx.constBegin(); it != ang_to_idx.constEnd(); ++it)
            {
                idx_to_idx.insert(it.value(), all_ang_to_idx.value(it.key(), -1));
            }

            // now run through the angles updating their local indicies to
            // global indicies
            const auto angs_inc_h = std::get<0>(all_ang_data_data[i]);
            const auto angs_exc_h = std::get<1>(all_ang_data_data[i]);

            for (int j = 0; j < angs_inc_h.count(); j += 4)
            {
                all_angs_inc_h.append(angs_inc_h[j]);
                all_angs_inc_h.append(angs_inc_h[j + 1]);
                all_angs_inc_h.append(angs_inc_h[j + 2]);
                all_angs_inc_h.append(idx_to_idx.value(angs_inc_h[j + 3]));
            }

            for (int j = 0; j < angs_exc_h.count(); j += 4)
            {
                all_angs_exc_h.append(angs_exc_h[j]);
                all_angs_exc_h.append(angs_exc_h[j + 1]);
                all_angs_exc_h.append(angs_exc_h[j + 2]);
                all_angs_exc_h.append(idx_to_idx.value(angs_exc_h[j + 3]));
            }
        }

        int nangh = all_angs_inc_h.count() / 4; // number of angles containing hydrogen
        int manga = all_angs_exc_h.count() / 4; // number of angles not containing hydrogen
        int nanga = manga;                      // manga + number of constraint angles
        int numang = k_data.count();            // number of unique angle types

        return std::make_tuple(writeFloatData(k_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeFloatData(t0_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeIntData(all_angs_inc_h, AmberFormat(AmberPrm::INTEGER, 10, 8)),
                               writeIntData(all_angs_exc_h, AmberFormat(AmberPrm::INTEGER, 10, 8)), nangh, manga, nanga,
                               numang);
    };

    // function used to get all of the dihedral information
    auto getAllDihedrals = [&]()
    {
        QVector<std::tuple<QVector<qint64>, QVector<qint64>, QHash<AmberNBDihPart, qint64>>> all_dih_data(
            params.count());

        auto all_dih_data_data = all_dih_data.data();

        // get all of the dihedral information from each molecule
        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    all_dih_data_data[i] = getDihedralData(params_data[i], atom_start_index_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                all_dih_data_data[i] = getDihedralData(params_data[i], atom_start_index_data[i]);
            }
        }

        // count up the number of dihedral parameters to help reserve memory space
        int ndihs_inc_h = 0;
        int ndihs_exc_h = 0;
        int nparams = 0;
        for (int i = 0; i < params.count(); ++i)
        {
            ndihs_inc_h += std::get<0>(all_dih_data_data[i]).count();
            ndihs_exc_h += std::get<1>(all_dih_data_data[i]).count();
            nparams += std::get<2>(all_dih_data_data[i]).count();
        }

        QVector<qint64> all_dihs_inc_h, all_dihs_exc_h;
        all_dihs_inc_h.reserve(ndihs_inc_h);
        all_dihs_exc_h.reserve(ndihs_exc_h);

        // we now need to go through the dihedrals and remove duplicated parameters
        QVector<double> force_data, per_data, phase_data, scee_data, scnb_data;
        QHash<AmberNBDihPart, int> all_dih_to_idx;
        all_dih_to_idx.reserve(nparams);

        // combine all parameters into a single set, with a unique index
        for (int i = 0; i < params.count(); ++i)
        {
            const auto dih_to_idx = std::get<2>(all_dih_data_data[i]);

            // first, find all unique dihedrals
            for (auto it = dih_to_idx.constBegin(); it != dih_to_idx.constEnd(); ++it)
            {
                if (not all_dih_to_idx.contains(it.key()))
                {
                    all_dih_to_idx.insert(it.key(), -1);
                }
            }
        }

        // now, sort the list so that there is a consistent order of
        // parameters every time this output file is written
        {
            auto all_dihs = all_dih_to_idx.keys();
            std::sort(all_dihs.begin(), all_dihs.end());

            force_data = QVector<double>(all_dihs.count());
            per_data = QVector<double>(all_dihs.count());
            phase_data = QVector<double>(all_dihs.count());
            scee_data = QVector<double>(all_dihs.count());
            scnb_data = QVector<double>(all_dihs.count());

            int i = 0;
            for (auto dih : all_dihs)
            {
                force_data[i] = dih.parameter().k();
                per_data[i] = dih.parameter().periodicity();
                phase_data[i] = dih.parameter().phase();

                if (dih.scl().cscl() != 0)
                {
                    scee_data[i] = 1.0 / dih.scl().cscl();
                }
                else
                {
                    scee_data[i] = 0.0;
                }

                if (dih.scl().ljscl() != 0)
                {
                    scnb_data[i] = 1.0 / dih.scl().ljscl();
                }
                else
                {
                    scnb_data[i] = 0.0;
                }

                all_dih_to_idx[dih] = i + 1;
                i += 1;
            }
        }

        for (int i = 0; i < params.count(); ++i)
        {
            // create a hash to map local parameter indicies to global parameter indicies
            QHash<qint64, qint64> idx_to_idx;

            // get all of the dihedral parameters
            const auto dih_to_idx = std::get<2>(all_dih_data_data[i]);

            // for each one, get the global parameter index, save the parameters to
            // the parameter arrays, and update the mapping from local to global indicies
            for (auto it = dih_to_idx.constBegin(); it != dih_to_idx.constEnd(); ++it)
            {
                idx_to_idx.insert(it.value(), all_dih_to_idx.value(it.key(), -1));
            }

            // now run through the dihedrals updating their local indicies to
            // global indicies
            const auto dihs_inc_h = std::get<0>(all_dih_data_data[i]);
            const auto dihs_exc_h = std::get<1>(all_dih_data_data[i]);

            for (int j = 0; j < dihs_inc_h.count(); j += 5)
            {
                all_dihs_inc_h.append(dihs_inc_h[j]);
                all_dihs_inc_h.append(dihs_inc_h[j + 1]);
                all_dihs_inc_h.append(dihs_inc_h[j + 2]);
                all_dihs_inc_h.append(dihs_inc_h[j + 3]);
                all_dihs_inc_h.append(idx_to_idx.value(dihs_inc_h[j + 4]));
            }

            for (int j = 0; j < dihs_exc_h.count(); j += 5)
            {
                all_dihs_exc_h.append(dihs_exc_h[j]);
                all_dihs_exc_h.append(dihs_exc_h[j + 1]);
                all_dihs_exc_h.append(dihs_exc_h[j + 2]);
                all_dihs_exc_h.append(dihs_exc_h[j + 3]);
                all_dihs_exc_h.append(idx_to_idx.value(dihs_exc_h[j + 4]));
            }
        }

        // now sort the arrays so that the order is the same every time this
        // molecule is written to a file
        ::detail::Idx5 *start_it = reinterpret_cast<::detail::Idx5 *>(all_dihs_inc_h.data());
        ::detail::Idx5 *end_it = start_it + (all_dihs_inc_h.count() / 5);
        std::sort(start_it, end_it);

        start_it = reinterpret_cast<::detail::Idx5 *>(all_dihs_exc_h.data());
        end_it = start_it + (all_dihs_exc_h.count() / 5);
        std::sort(start_it, end_it);

        auto ignore_14 = [](QVector<qint64> &dihs)
        {
            QSet<QPair<qint64, qint64>> seen_dihedrals;
            seen_dihedrals.reserve(dihs.count() / 5);

            for (int i = 0; i < dihs.count(); i += 5)
            {
                // NB14 must be computed only the first dihedral found that has this pair,
                // otherwise we would risk Amber double-counting this parameter
                QPair<qint64, qint64> nb14pair(dihs[i], dihs[i + 3]);

                if (dihs[i] > dihs[i + 3])
                {
                    // must ensure a consistent ordering, so that we catch
                    // A-*-*-B and B-*-*-A
                    nb14pair = QPair<qint64, qint64>(dihs[i + 3], dihs[i]);
                }

                if (seen_dihedrals.contains(nb14pair))
                {
                    dihs[i + 2] = -dihs[i + 2];
                }
                else if (dihs[i + 2] > 0)
                {
                    seen_dihedrals.insert(nb14pair);
                }
            }
        };

        ignore_14(all_dihs_inc_h);
        ignore_14(all_dihs_exc_h);

        int ndihh = all_dihs_inc_h.count() / 5; // number of dihedrals containing hydrogen
        int mdiha = all_dihs_exc_h.count() / 5; // number of dihedrals not containing hydrogen
        int ndiha = mdiha;                      // mdiha + number of constraint dihedrals
        int numdih = force_data.count();        // number of unique dihedral types

        return std::make_tuple(writeFloatData(force_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeFloatData(per_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeFloatData(phase_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeFloatData(scee_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeFloatData(scnb_data, AmberFormat(AmberPrm::FLOAT, 5, 16, 8)),
                               writeIntData(all_dihs_inc_h, AmberFormat(AmberPrm::INTEGER, 10, 8)),
                               writeIntData(all_dihs_exc_h, AmberFormat(AmberPrm::INTEGER, 10, 8)), ndihh, mdiha, ndiha,
                               numdih);
    };

    // function used to get all of the CMAP information
    auto getAllCMAPs = [&]()
    {
        QVector<std::tuple<QVector<qint64>, QHash<CMAPParameter, qint64>>> all_cmap_data(params.count());
        auto all_cmap_data_data = all_cmap_data.data();

        // get all of the cmap information from each molecule
        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, params.count()), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    all_cmap_data_data[i] = getCMAPData(params_data[i], atom_start_index_data[i]);
                } });
        }
        else
        {
            for (int i = 0; i < params.count(); ++i)
            {
                all_cmap_data_data[i] = getCMAPData(params_data[i], atom_start_index_data[i]);
            }
        }

        // count up the number of cmap parameters to help reserve memory space
        int ncmaps = 0;
        int nparams = 0;
        for (int i = 0; i < params.count(); ++i)
        {
            ncmaps += std::get<0>(all_cmap_data[i]).count();
            nparams += std::get<1>(all_cmap_data[i]).count();
        }

        QVector<qint64> all_cmap_idxs;
        all_cmap_idxs.reserve(ncmaps);

        // we now need to go through the cmaps and remove duplicated cmap parameters
        QVector<CMAPParameter> cmap_params;
        QHash<CMAPParameter, int> all_cmap_to_idx;
        all_cmap_to_idx.reserve(nparams);

        // combine all parameters into a single set, with a unique index
        for (int i = 0; i < params.count(); ++i)
        {
            const auto cmap_to_idx = std::get<1>(all_cmap_data_data[i]);

            // first, find all unique cmaps
            for (auto it = cmap_to_idx.constBegin(); it != cmap_to_idx.constEnd(); ++it)
            {
                if (not all_cmap_to_idx.contains(it.key()))
                {
                    all_cmap_to_idx.insert(it.key(), -1);
                }
            }
        }

        // now, sort the list so that there is a consistent order of
        // parameters every time this output file is written
        {
            auto all_cmaps = all_cmap_to_idx.keys();
            std::sort(all_cmaps.begin(), all_cmaps.end());

            cmap_params = QVector<CMAPParameter>(all_cmaps.count());

            int i = 0;
            for (const auto &cmap : all_cmaps)
            {
                cmap_params[i] = cmap;

                all_cmap_to_idx[cmap] = i + 1;
                i += 1;
            }
        }

        for (int i = 0; i < params.count(); ++i)
        {
            // create a hash to map local parameter indicies to global parameter indicies
            QHash<qint64, qint64> idx_to_idx;

            // get all of the cmap parameters
            const auto cmap_to_idx = std::get<1>(all_cmap_data_data[i]);

            // for each one, get the global parameter index, save the parameters to
            // the cmap arrays, and update the mapping from local to global indicies
            for (auto it = cmap_to_idx.constBegin(); it != cmap_to_idx.constEnd(); ++it)
            {
                idx_to_idx.insert(it.value(), all_cmap_to_idx.value(it.key(), -1));
            }

            // now run through the cmaps updating their local indicies to
            // global indicies
            const auto cmap_idxs = std::get<0>(all_cmap_data_data[i]);

            for (int j = 0; j < cmap_idxs.count(); j += 6)
            {
                all_cmap_idxs.append(cmap_idxs[j]);
                all_cmap_idxs.append(cmap_idxs[j + 1]);
                all_cmap_idxs.append(cmap_idxs[j + 2]);
                all_cmap_idxs.append(cmap_idxs[j + 3]);
                all_cmap_idxs.append(cmap_idxs[j + 4]);
                all_cmap_idxs.append(idx_to_idx.value(cmap_idxs[j + 5]));
            }
        }

        int ncmap = all_cmap_idxs.count() / 6; // number of cmap terms

        return std::make_tuple(writeIntData(all_cmap_idxs, AmberFormat(AmberPrm::INTEGER, 6, 8)),
                               cmap_params, ncmap);
    };

    QStringList lines;

    QStringList name_lines, charge_lines, number_lines, mass_lines, radius_lines, ambtyp_lines, screening_lines,
        treechain_lines;
    std::tuple<QStringList, QStringList, QStringList, QStringList, QStringList> lj_lines;
    std::tuple<QStringList, QStringList, int> excl_lines;
    std::tuple<QStringList, QStringList, int, int> res_lines;
    std::tuple<QStringList, QStringList, QStringList, QStringList, int, int, int, int> bond_lines;
    std::tuple<QStringList, QStringList, QStringList, QStringList, int, int, int, int> ang_lines;
    std::tuple<QStringList, QStringList, QStringList, QStringList, QStringList, QStringList, QStringList, int, int, int,
               int>
        dih_lines;
    std::tuple<QStringList, QVector<CMAPParameter>, int> cmap_lines;

    const QVector<std::function<void()>> info_functions = {
        [&]()
        { excl_lines = getAllExcludedAtoms(); },
        [&]()
        { bond_lines = getAllBonds(); },
        [&]()
        { ang_lines = getAllAngles(); },
        [&]()
        { dih_lines = getAllDihedrals(); },
        [&]()
        { cmap_lines = getAllCMAPs(); },
        [&]()
        { name_lines = getAllAtomNames(); },
        [&]()
        { charge_lines = getAllAtomCharges(); },
        [&]()
        { number_lines = getAllAtomNumbers(); },
        [&]()
        { mass_lines = getAllAtomMasses(); },
        [&]()
        { lj_lines = getAllAtomTypes(); },
        [&]()
        { ambtyp_lines = getAllAmberTypes(); },
        [&]()
        { radius_lines = getAllRadii(); },
        [&]()
        { screening_lines = getAllScreening(); },
        [&]()
        { treechain_lines = getAllTreeChains(); },
        [&]()
        { res_lines = getAllResidueInfo(); }};

    if (use_parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, info_functions.count()), [&](const tbb::blocked_range<int> &r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                info_functions[i]();
            } });
    }
    else
    {
        for (auto info_function : info_functions)
        {
            info_function();
        }
    }

    // save the number of excluded atoms to 'pointers'
    pointers[10] = std::get<2>(excl_lines); // NNB

    // save the number of residues to 'pointers'
    pointers[11] = std::get<2>(res_lines); // NRES
    pointers[28] = std::get<3>(res_lines); // NMXRS

    // save the number of bonds to 'pointers'
    pointers[2] = std::get<4>(bond_lines);  // NBONH
    pointers[3] = std::get<5>(bond_lines);  // MBONA
    pointers[12] = std::get<6>(bond_lines); // NBONA
    pointers[15] = std::get<7>(bond_lines); // NUMBND

    // save the number of angles to 'pointers'
    pointers[4] = std::get<4>(ang_lines);  // NTHETH
    pointers[5] = std::get<5>(ang_lines);  // MTHETA
    pointers[13] = std::get<6>(ang_lines); // NTHETA
    pointers[16] = std::get<7>(ang_lines); // NUMANG

    // save the number of dihedrals to 'pointers'
    pointers[6] = std::get<7>(dih_lines);   // NPHIH
    pointers[7] = std::get<8>(dih_lines);   // MPHIA
    pointers[14] = std::get<9>(dih_lines);  // NPHIA
    pointers[17] = std::get<10>(dih_lines); // NPTRA

    // Add the number of dummy atoms, i.e. NUMEXTRA.
    pointers[30] = num_dummies;

    const int natyp = pointers[18]; // NATYP

    lines.append("%FLAG POINTERS");
    lines += writeIntData(pointers, AmberFormat(AmberPrm::INTEGER, 10, 8));

    lines.append("%FLAG ATOM_NAME");
    lines += name_lines;

    lines.append("%FLAG CHARGE");
    lines += charge_lines;

    lines.append("%FLAG ATOMIC_NUMBER");
    lines += number_lines;

    lines.append("%FLAG MASS");
    lines += mass_lines;

    lines.append("%FLAG ATOM_TYPE_INDEX");
    lines += std::get<0>(lj_lines);

    lines.append("%FLAG NUMBER_EXCLUDED_ATOMS");
    lines += std::get<0>(excl_lines);

    lines.append("%FLAG NONBONDED_PARM_INDEX");
    lines += std::get<1>(lj_lines);

    lines.append("%FLAG RESIDUE_LABEL");
    lines += std::get<0>(res_lines);

    lines.append("%FLAG RESIDUE_POINTER");
    lines += std::get<1>(res_lines);

    lines.append("%FLAG BOND_FORCE_CONSTANT");
    lines += std::get<0>(bond_lines);

    lines.append("%FLAG BOND_EQUIL_VALUE");
    lines += std::get<1>(bond_lines);

    lines.append("%FLAG ANGLE_FORCE_CONSTANT");
    lines += std::get<0>(ang_lines);

    lines.append("%FLAG ANGLE_EQUIL_VALUE");
    lines += std::get<1>(ang_lines);

    lines.append("%FLAG DIHEDRAL_FORCE_CONSTANT");
    lines += std::get<0>(dih_lines);

    lines.append("%FLAG DIHEDRAL_PERIODICITY");
    lines += std::get<1>(dih_lines);

    lines.append("%FLAG DIHEDRAL_PHASE");
    lines += std::get<2>(dih_lines);

    lines.append("%FLAG SCEE_SCALE_FACTOR");
    lines += std::get<3>(dih_lines);

    lines.append("%FLAG SCNB_SCALE_FACTOR");
    lines += std::get<4>(dih_lines);

    lines.append("%FLAG SOLTY");
    // this is currently unused in Amber and should be equal to 0.0 for all atom types
    {
        QVector<double> solty(natyp, 0.0);
        lines += writeFloatData(solty, AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
    }

    lines.append("%FLAG LENNARD_JONES_ACOEF");
    lines += std::get<2>(lj_lines);

    lines.append("%FLAG LENNARD_JONES_BCOEF");
    lines += std::get<3>(lj_lines);

    if (not std::get<4>(lj_lines).isEmpty())
    {
        lines.append("%FLAG LENNARD_JONES_CCOEF");
        lines += std::get<4>(lj_lines);
    }

    lines.append("%FLAG BONDS_INC_HYDROGEN");
    lines += std::get<2>(bond_lines);

    lines.append("%FLAG BONDS_WITHOUT_HYDROGEN");
    lines += std::get<3>(bond_lines);

    lines.append("%FLAG ANGLES_INC_HYDROGEN");
    lines += std::get<2>(ang_lines);

    lines.append("%FLAG ANGLES_WITHOUT_HYDROGEN");
    lines += std::get<3>(ang_lines);

    lines.append("%FLAG DIHEDRALS_INC_HYDROGEN");
    lines += std::get<5>(dih_lines);

    lines.append("%FLAG DIHEDRALS_WITHOUT_HYDROGEN");
    lines += std::get<6>(dih_lines);

    lines.append("%FLAG EXCLUDED_ATOMS_LIST");
    lines += std::get<1>(excl_lines);

    // HBOND_ACOEF parameters are not currently supported
    lines.append("%FLAG HBOND_ACOEF");
    lines.append("%FORMAT(5E16.8)");
    lines.append(" ");

    // HBOND_BCOEF parameters are not currently supported
    lines.append("%FLAG HBOND_BCOEF");
    lines.append("%FORMAT(5E16.8)");
    lines.append(" ");

    // HBCUT parameters not supported, but now not used in Amber
    lines.append("%FLAG HBCUT");
    lines.append("%FORMAT(5E16.8)");
    lines.append(" ");

    lines.append("%FLAG AMBER_ATOM_TYPE");
    lines += ambtyp_lines;

    lines.append("%FLAG TREE_CHAIN_CLASSIFICATION");
    lines += treechain_lines;

    lines.append("%FLAG JOIN_ARRAY");
    {
        // no longer used in Amber, and should be an array of zeroes
        QVector<qint64> join_array(total_natoms, 0);
        auto zeroes = writeIntData(join_array, AmberFormat(AmberPrm::INTEGER, 10, 8));
        lines += zeroes;

        // IROTAT is also no longer used, and is also an array of zeroes
        lines.append("%FLAG IROTAT");
        lines += zeroes;
    }

    lines.append("%FLAG RADIUS_SET");
    lines.append("%FORMAT(1a80)");
    // have to take the first molecule's radius set as the PRM file assumes that
    // they are all the same...
    lines.append(params[0].radiusSet());

    lines.append("%FLAG RADII");
    lines += radius_lines;

    lines.append("%FLAG SCREEN");
    lines += screening_lines;

    lines.append("%FLAG ATOMS_PER_MOLECULE");
    lines += writeIntData(natoms_per_molecule, AmberFormat(AmberPrm::INTEGER, 10, 8));

    if (has_periodic_box)
    {
        // write out the "SOLVENT_POINTERS". Amber has a concept of solute(s) and solvent.
        // with all solutes written first, and then all solvents. We have above calculated
        // the number of solvents by finding the last "nsolvents" molecules that have the
        // same number of atoms, with the same atom names in the same order. We are relying
        // on being passed a pre-sorted system that has all solvent molecules at the end
        lines.append("%FLAG SOLVENT_POINTERS");
        QVector<qint64> solvent_pointers(3, 0);

        if (nsolvents > 0)
        {
            solvent_pointers[0] = last_solute_residue;            // last solute residue index
            solvent_pointers[1] = params.count();                 // total number of molecules
            solvent_pointers[2] = params.count() - nsolvents + 1; // first solvent molecule index
        }

        lines += writeIntData(solvent_pointers, AmberFormat(AmberPrm::INTEGER, 3, 8));

        // write out the box dimensions
        if (space.isA<PeriodicBox>())
        {
            const auto dims = space.asA<PeriodicBox>().dimensions();
            QVector<double> box_dims(4);

            box_dims[0] = 90.0;
            box_dims[1] = dims.x();
            box_dims[2] = dims.y();
            box_dims[3] = dims.z();

            lines.append("%FLAG BOX_DIMENSIONS");
            lines += writeFloatData(box_dims, AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
        }
        else if (space.isA<TriclinicBox>())
        {
            const auto v0 = space.asA<TriclinicBox>().vector0();
            const auto v1 = space.asA<TriclinicBox>().vector1();
            const auto v2 = space.asA<TriclinicBox>().vector2();
            double beta = space.asA<TriclinicBox>().beta();

            QVector<double> box_dims(4);

            QVector<double> boxdims(6);
            box_dims[0] = beta;
            box_dims[1] = v0.magnitude();
            box_dims[2] = v1.magnitude();
            box_dims[3] = v2.magnitude();

            lines.append("%FLAG BOX_DIMENSIONS");
            lines += writeFloatData(box_dims, AmberFormat(AmberPrm::FLOAT, 5, 16, 8));
        }
    }

    // we currently only support fixed-charge forcefields (IPOL = 0)
    lines.append("%FLAG IPOL");
    lines.append("%FORMAT(1I8)");
    lines.append("       0");

    // we don't currently support IFCAP > 0, IFPERT > 0 or IFPOL > 0

    if (std::get<2>(cmap_lines) > 0)
    {
        const auto cmap_params = std::get<1>(cmap_lines);

        lines.append("%FLAG CMAP_COUNT");
        QVector<qint64> counts(2);
        counts[0] = std::get<2>(cmap_lines);
        counts[1] = cmap_params.count();
        lines += writeIntData(counts, AmberFormat(AmberPrm::INTEGER, 2, 8));

        QVector<qint64> resolutions;

        for (const auto &cmap_param : cmap_params)
        {
            if (cmap_param.nRows() != cmap_param.nColumns())
            {
                throw SireIO::parse_error(
                    QObject::tr("CMAP parameters must be square matrices, but the matrix for %1x%2")
                        .arg(cmap_param.nRows())
                        .arg(cmap_param.nColumns()),
                    CODELOC);
            }

            resolutions.append(cmap_param.nRows());
        }

        lines.append("%FLAG CMAP_RESOLUTION");
        lines += writeIntData(resolutions, AmberFormat(AmberPrm::INTEGER, 20, 4));

        for (int i = 0; i < cmap_params.count(); ++i)
        {
            const auto &cmap_param = cmap_params[i];

            lines.append(QString("%FLAG CMAP_PARAMETER_%1").arg(i + 1, 2, 10, QChar('0')));

            lines += writeFloatData(cmap_param.grid().toColumnMajorVector(),
                                    AmberFormat(AmberPrm::FFLOAT, 8, 9, 5));
        }

        lines.append("%FLAG CMAP_INDEX");
        lines += std::get<0>(cmap_lines);
    }

    return lines;
}

/** Construct by converting from the passed system, using the passed property
    map to find the right properties */
AmberPrm::AmberPrm(const System &system, const PropertyMap &map) : ConcreteProperty<AmberPrm, MoleculeParser>(map)
{
    // get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();
    const auto molnums_data = molnums.constData();

    if (molnums.isEmpty())
    {
        // no molecules in the system
        this->operator=(AmberPrm());
        return;
    }

    try
    {
        ffield = system.property(map["forcefield"]).asA<MMDetail>();
    }
    catch (...)
    {
        ffield = MMDetail::guessFrom("arithmetic", "coulomb", "lj", 1.0 / 1.2, 0.5, "harmonic", "harmonic", "cosine");
    }

    if (ffield.usesGeometricCombiningRules())
    {
        comb_rules = "geometric";
    }
    else
    {
        comb_rules = "arithmetic";
    }

    // generate the AmberParams object for each molecule in the system
    QVector<AmberParams> params(molnums.count());
    auto params_data = params.data();

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, molnums.count()), [&](const tbb::blocked_range<int> r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                params_data[i] = AmberParams(system[molnums_data[i]], map);
            } });
    }
    else
    {
        for (int i = 0; i < molnums.count(); ++i)
        {
            params_data[i] = AmberParams(system[molnums_data[i]], map);
        }
    }

    QStringList errors;

    // extract the space of the system
    SpacePtr space;

    try
    {
        space = system.property(map["space"]).asA<Space>();
    }
    catch (...)
    {
    }

    // Work out the number of dummy atoms to add the NUMEXTRA pointer.
    int num_dummies = system.search("element Xx").count();

    // now convert these into text lines that can be written as the file
    QStringList lines = ::toLines(params, space, num_dummies,
                                  this->comb_rules == "geometric",
                                  this->usesParallel(), &errors);

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error(
            QObject::tr("Errors converting the system to a Amber Parm format...\n%1").arg(errors.join("\n")), CODELOC);
    }

    // we don't need params any more, so free the memory
    params.clear();

    // add the title and the date to the top of the lines
    lines.prepend(system.name().value());
    lines.prepend("%FORMAT(20a4)");
    lines.prepend("%FLAG TITLE");
    lines.prepend(QString("%VERSION  VERSION_STAMP = V0001.000  DATE = %1")
                      .arg(QDateTime::currentDateTime().toString("MM/dd/yy  hh:mm:ss")));

    // now generate this object by re-reading these lines
    AmberPrm parsed(lines, map);

    this->operator=(parsed);
}

/** Copy constructor */
AmberPrm::AmberPrm(const AmberPrm &other)
    : ConcreteProperty<AmberPrm, MoleculeParser>(other), flag_to_line(other.flag_to_line), int_data(other.int_data),
      float_data(other.float_data), string_data(other.string_data), lj_data(other.lj_data),
      lj_exceptions(other.lj_exceptions), cmap_data(other.cmap_data),
      bonds_inc_h(other.bonds_inc_h), bonds_exc_h(other.bonds_exc_h), angs_inc_h(other.angs_inc_h),
      angs_exc_h(other.angs_exc_h), dihs_inc_h(other.dihs_inc_h), dihs_exc_h(other.dihs_exc_h),
      excl_atoms(other.excl_atoms), cmap_idxs(other.cmap_idxs),
      molnum_to_atomnums(other.molnum_to_atomnums), pointers(other.pointers),
      ffield(other.ffield), warns(other.warns), comb_rules(other.comb_rules)
{
}

/** Destructor */
AmberPrm::~AmberPrm()
{
}

/** Copy assignment operator */
AmberPrm &AmberPrm::operator=(const AmberPrm &other)
{
    if (this != &other)
    {
        flag_to_line = other.flag_to_line;
        int_data = other.int_data;
        float_data = other.float_data;
        string_data = other.string_data;
        lj_data = other.lj_data;
        lj_exceptions = other.lj_exceptions;
        cmap_data = other.cmap_data;
        bonds_inc_h = other.bonds_inc_h;
        bonds_exc_h = other.bonds_exc_h;
        angs_inc_h = other.angs_inc_h;
        angs_exc_h = other.angs_exc_h;
        dihs_inc_h = other.dihs_inc_h;
        dihs_exc_h = other.dihs_exc_h;
        excl_atoms = other.excl_atoms;
        cmap_idxs = other.cmap_idxs;
        molnum_to_atomnums = other.molnum_to_atomnums;
        pointers = other.pointers;
        ffield = other.ffield;
        warns = other.warns;
        comb_rules = other.comb_rules;
        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool AmberPrm::operator==(const AmberPrm &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool AmberPrm::operator!=(const AmberPrm &other) const
{
    return not operator==(other);
}

const char *AmberPrm::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberPrm>());
}

const char *AmberPrm::what() const
{
    return AmberPrm::typeName();
}

/** Return a string representation of this object */
QString AmberPrm::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("AmberPrm::null");
    else
    {
        return QObject::tr("AmberPrm( title() = '%1', nMolecules() = %2, "
                           "nResidues() = %3, nAtoms() = %4 )")
            .arg(title())
            .arg(nMolecules())
            .arg(nResidues())
            .arg(nAtoms());
    }
}

QStringList AmberPrm::warnings() const
{
    return warns;
}

/** Return the title of the parameter file */
QString AmberPrm::title() const
{
    return QStringList(string_data.value("TITLE").toList()).join("").simplified();
}

/** Return the total number of atoms in the file */
int AmberPrm::nAtoms() const
{
    if (pointers.count() > 0)
        return pointers[0];
    else
        return 0;
}

/** Return the total number of atoms in the ith molecule in the file */
int AmberPrm::nAtoms(int idx) const
{
    idx = Index(idx).map(this->nMolecules());
    return molnum_to_atomnums.at(idx + 1).count(); // molnum is 1-indexed
}

/** Return the number of distinct atom types */
int AmberPrm::nTypes() const
{
    if (pointers.count() > 1)
        return pointers[1];
    else
        return 0;
}

/** Return the total number of bonds */
int AmberPrm::nBonds() const
{
    return nBondsWithHydrogen() + nBondsNoHydrogen();
}

/** Return the number of bonds including hydrogen */
int AmberPrm::nBondsWithHydrogen() const
{
    if (pointers.count() > 2)
        return pointers[2];
    else
        return 0;
}

/** Return the number of bonds no containing hydrogen */
int AmberPrm::nBondsNoHydrogen() const
{
    if (pointers.count() > 3)
        return pointers[3];
    else
        return 0;
}

/** Return the number of angles */
int AmberPrm::nAngles() const
{
    return nAnglesWithHydrogen() + nAnglesNoHydrogen();
}

/** Return the number of angles containing hydrogen */
int AmberPrm::nAnglesWithHydrogen() const
{
    if (pointers.count() > 4)
        return pointers[4];
    else
        return 0;
}

/** Return the number of angles without hydrogen */
int AmberPrm::nAnglesNoHydrogen() const
{
    if (pointers.count() > 5)
        return pointers[5];
    else
        return 0;
}

/** Return the number of dihedrals */
int AmberPrm::nDihedrals() const
{
    return nDihedralsWithHydrogen() + nDihedralsNoHydrogen();
}

/** Return the number of dihedrals containing hydrogen */
int AmberPrm::nDihedralsWithHydrogen() const
{
    if (pointers.count() > 6)
        return pointers[6];
    else
        return 0;
}

/** Return the number of dihedrals without hydrogen */
int AmberPrm::nDihedralsNoHydrogen() const
{
    if (pointers.count() > 7)
        return pointers[7];
    else
        return 0;
}

/** Return the number of excluded atoms */
int AmberPrm::nExcluded() const
{
    if (pointers.count() > 10)
        return pointers[10];
    else
        return 0;
}

/** Return the number of residues */
int AmberPrm::nResidues() const
{
    if (pointers.count() > 11)
        return pointers[11];
    else
        return 0;
}

/** Return the number of molecules in the file */
int AmberPrm::nMolecules() const
{
    QVector<qint64> atoms_per_mol = int_data.value("ATOMS_PER_MOLECULE");

    if (atoms_per_mol.isEmpty())
    {
        if (nAtoms() == 0)
            return 0;
        else
            return 1;
    }
    else
        return atoms_per_mol.count();
}

/** Return an array that maps from the index of each atom to the index
    of each molecule. Note that both the atom index and molecule index
    are 0-indexed */
QVector<int> AmberPrm::getAtomIndexToMolIndex() const
{
    const int natoms = this->nAtoms();

    if (natoms <= 0)
        return QVector<int>();

    QVector<int> atom_to_mol(natoms);
    auto atom_to_mol_array = atom_to_mol.data();

    const int nmols = molnum_to_atomnums.count() - 1; // 1 indexed

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                const auto atoms_in_mol = molnum_to_atomnums.at(i + 1);

                for (const auto atom : atoms_in_mol)
                {
                    atom_to_mol_array[atom - 1] = i;
                }
            } });
    }
    else
    {
        for (int i = 0; i < nmols; ++i)
        {
            const auto atoms_in_mol = molnum_to_atomnums.at(i + 1);

            for (const auto atom : atoms_in_mol)
            {
                atom_to_mol_array[atom - 1] = i;
            }
        }
    }

    return atom_to_mol;
}

/** Return an AmberPrm object parsed from the passed file */
AmberPrm AmberPrm::parse(const QString &filename, const PropertyMap &map)
{
    return AmberPrm(filename);
}

/** Return the lines that correspond to the passed flag. This returns an
    empty list of there are no lines associated with the passed flag */
QVector<QString> AmberPrm::linesForFlag(const QString &flag) const
{
    auto it = flag_to_line.constFind(flag);

    if (it != flag_to_line.constEnd())
    {
        const int start = it->first;
        const int count = it->second;

        SireBase::assert_true(start >= 0 and start < lines().count(), CODELOC);
        SireBase::assert_true(count > 0 and start + count < lines().count(), CODELOC);

        return lines().mid(start, count);
    }
    else
        return QVector<QString>();
}

/** Return all of the flags that are held in this file */
QStringList AmberPrm::flags() const
{
    return flag_to_line.keys();
}

template <class FUNC, class Exception>
void live_test(FUNC function, QList<std::shared_ptr<Exception>> &errors)
{
    try
    {
        return function();
    }
    catch (const Exception &e)
    {
        errors.append(std::shared_ptr<Exception>(e.clone()));
    }
}

/** Run through all of the data that has been read and perform a series
    of tests that will see if the prm7 data is sane. If any test fails,
    then an exception will be thrown */
void AmberPrm::assertSane() const
{
    QList<std::shared_ptr<SireError::exception>> errors;

    int natoms = this->nAtoms();

    live_test([&]()
              { SireBase::assert_equal(this->stringData("ATOM_NAME").count(), natoms, CODELOC); },
              errors);

    live_test([&]()
              { SireBase::assert_equal(this->floatData("CHARGE").count(), natoms, CODELOC); },
              errors);

    live_test([&]()
              { SireBase::assert_equal(this->intData("ATOMIC_NUMBER").count(), natoms, CODELOC); },
              errors);

    live_test([&]()
              { SireBase::assert_equal(this->floatData("MASS").count(), natoms, CODELOC); },
              errors);

    live_test([&]()
              { SireBase::assert_equal(this->intData("ATOM_TYPE_INDEX").count(), natoms, CODELOC); },
              errors);

    live_test([&]()
              { SireBase::assert_equal(this->floatData("RADII").count(), natoms, CODELOC); },
              errors);

    live_test([&]()
              { SireBase::assert_equal(this->floatData("SCREEN").count(), natoms, CODELOC); },
              errors);

    live_test([&]()
              { SireBase::assert_equal(this->stringData("TREE_CHAIN_CLASSIFICATION").count(), natoms, CODELOC); },
              errors);

    if (not errors.isEmpty())
    {
        for (auto error : errors)
        {
            qDebug() << error->toString();
        }

        throw SireIO::parse_error(QObject::tr("Sanity tests failed for the loaded Amber prm7 format file"), CODELOC);
    }
}

/** Internal function used to get the molecule structure that starts at index 'start_idx'
    in the file, and that has 'natoms' atoms, and is molecule at index 'molidx'. Note
    that if this 'molecule' is actually bonded to another 'molecule' in the input file,
    then we will only return a MolStructureEditor for the combined molecule with the
    lowest molidx - a null MolStructureEditor will be returned for the other molecules
*/
MolEditor AmberPrm::getMolStructure(int molidx, const PropertyName &cutting) const
{
    const int molnum = molidx + 1; // amber is 1-indexed

    const auto atomnums_in_mol = molnum_to_atomnums.at(molnum);

    const int natoms = atomnums_in_mol.count();

    if (natoms == 0)
    {
        throw SireError::program_bug(QObject::tr("Strange - there are no atoms in this molecule %1?").arg(molnum),
                                     CODELOC);
    }

    // locate the residue pointers for this molecule - note that the
    // residue pointers are 1-indexed (i.e. from atom 1 to atom N)
    const auto res_pointers = this->intData("RESIDUE_POINTER");

    QSet<int> res_idxs_set;

    int highest_atomnum = atomnums_in_mol[0];
    int lowest_atomnum = highest_atomnum;

    for (int i = 0; i < atomnums_in_mol.count(); ++i)
    {
        auto atomnum = atomnums_in_mol[i];

        if (atomnum > highest_atomnum)
        {
            highest_atomnum = atomnum;
        }

        if (atomnum < lowest_atomnum)
        {
            lowest_atomnum = atomnum;
        }
    }

    for (int i = 0; i < res_pointers.count(); ++i)
    {
        const int res_start_atom = res_pointers[i];

        if (res_start_atom >= lowest_atomnum)
        {
            if (atomnums_in_mol.contains(res_start_atom))
            {
                res_idxs_set.insert(i);
            }
        }
        else if (res_start_atom > highest_atomnum)
        {
            break;
        }
    }

    QList<int> res_idxs = QList<int>(res_idxs_set.constBegin(),
                                     res_idxs_set.constEnd());
    std::sort(res_idxs.begin(), res_idxs.end());

    if (res_idxs.isEmpty())
    {
        // no residues?
        throw SireError::program_bug(QObject::tr("Strange - there are no residues in this molecule %1?").arg(molnum),
                                     CODELOC);
    }

    auto raw_atom_names = this->stringData("ATOM_NAME");

    while (raw_atom_names.count() < highest_atomnum)
    {
        // we are missing some atom names - give a default
        // name of "UNK", which tends to be used for "UNKNOWN" atoms.
        // see
        // https://www.rcsb.org/ligand/UNX
        // for this convention
        raw_atom_names.append("UNK");
    }

    // make sure we use a const reference, or else we will
    // severely degrade performance by continually copying
    // this list of names
    const auto atom_names = raw_atom_names;

    auto raw_res_names = this->stringData("RESIDUE_LABEL");

    while (raw_res_names.count() <= res_idxs.at(res_idxs.count() - 1))
    {
        // we are missing some residue names - give a default
        // name of "UNK", which tends to be used for "UNKNOWN" residues
        // see
        // http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_37.html
        // for this convention
        raw_res_names.append("UNK");
    }

    const auto res_names = raw_res_names;

    // do this a quicker way if this is a single-residue molecule
    // (e.g. like water)
    if (res_idxs.count() == 1)
    {
        int res_idx = res_idxs[0];

        int res_start_atom = res_pointers[res_idx];
        int res_end_atom;

        if (res_idx + 1 == res_pointers.count())
        {
            res_end_atom = this->nAtoms(); // 1-indexed
        }
        else
        {
            res_end_atom = res_pointers[res_idx + 1] - 1; // 1 lower than first index of atom in next residue
        }

        qint64 resnum = res_idx + 1;
        QString resname = res_names[res_idx];

        const int nres_atoms = res_end_atom - res_start_atom + 1;

        QVector<QString> res_atomnames(nres_atoms);
        QVector<qint64> res_atomnums(nres_atoms);

        auto res_atomnames_data = res_atomnames.data();
        auto res_atomnums_data = res_atomnums.data();

        int atomidx = 0;

        for (int i = res_start_atom; i <= res_end_atom; ++i)
        {
            res_atomnames_data[atomidx] = atom_names[i - 1].trimmed();
            res_atomnums_data[atomidx] = i;
            atomidx += 1;
        }

        auto mol = Molecule(resname, MoleculeInfo(resname, resnum,
                                                  res_atomnames,
                                                  res_atomnums));

        return mol.edit();
    }

    // first step is to build the structure of the molecule - i.e.
    // the layout of cutgroups, residues and atoms
    MolStructureEditor mol;

    for (int i = 0; i < res_idxs.count(); ++i)
    {
        int res_idx = res_idxs.at(i);
        int res_num = res_idx + 1;
        int res_start_atom = res_pointers[res_idx];

        if (atomnums_in_mol.contains(res_start_atom))
        {
            // find the number of the last atom in the residue (which should be in the molecule?)
            int res_end_atom;

            if (res_idx + 1 == res_pointers.count())
            {
                res_end_atom = this->nAtoms(); // 1-indexed
            }
            else
            {
                res_end_atom = res_pointers[res_idx + 1] - 1; // 1 lower than first index
                                                              // of atom in next residue
            }

            if (not atomnums_in_mol.contains(res_end_atom))
            {
                throw SireIO::parse_error(QObject::tr("Atoms in residue %1 are not listed as part of the expected "
                                                      "molecule %2. Should be atoms %3 to %4, but molecule has %5.")
                                              .arg(res_num)
                                              .arg(molnum)
                                              .arg(res_start_atom)
                                              .arg(res_end_atom)
                                              .arg(Sire::toString(atomnums_in_mol)),
                                          CODELOC);
            }

            // we assume that this molecules contains the entire residue
            auto res = mol.add(ResNum(res_num));
            res.rename(ResName(res_names[res_idx].trimmed()));

            // by default we will use one CutGroup per residue - this
            // may be changed later by the cutting system.
            auto cutgroup = mol.add(CGName(QString::number(i)));

            for (int j = res_start_atom; j <= res_end_atom; ++j)
            {
                if (not atomnums_in_mol.contains(j))
                {
                    throw SireIO::parse_error(
                        QObject::tr("Atom %6 in residue %1 are not listed as part of the expected "
                                    "molecule %2. Should be atoms %3 to %4, but molecule has %5.")
                            .arg(res_num)
                            .arg(molnum)
                            .arg(res_start_atom)
                            .arg(res_end_atom)
                            .arg(Sire::toString(atomnums_in_mol))
                            .arg(j),
                        CODELOC);
                }
                else
                {
                    auto atom = cutgroup.add(AtomNum(j));
                    atom.rename(AtomName(atom_names[j - 1].trimmed())); // 0-index versus 1-index
                    atom.reparent(ResNum(res_num));
                }
            }
        }
    }

    if (cutting.hasValue())
    {
        const CuttingFunction &cutfunc = cutting.value().asA<CuttingFunction>();

        if (not cutfunc.isA<ResidueCutting>())
        {
            mol = cutfunc(mol);
        }
    }

    return mol.commit().edit();
}

/** Return the AmberParams for the ith molecule */
AmberParams AmberPrm::getAmberParams(int molidx, const MoleculeInfoData &molinfo) const
{
    // function to assign all of the atom parameters
    auto assign_atoms = [&]()
    {
        // assign all of the atomic properties
        const auto charge_array = this->floatData("CHARGE").constData();
        const auto mass_array = this->floatData("MASS").constData();
        const auto atomic_num_array = this->intData("ATOMIC_NUMBER").constData();
        const auto amber_type_array = this->intData("ATOM_TYPE_INDEX").constData();
        const auto ambertype_array = this->stringData("AMBER_ATOM_TYPE").constData();
        const auto born_radii_array = this->floatData("RADII").constData();
        const auto born_screening_array = this->floatData("SCREEN").constData();
        const auto treechains_array = this->stringData("TREE_CHAIN_CLASSIFICATION").constData();

        const int natoms = molinfo.nAtoms();

        AmberParams params(molinfo);

        for (int i = 0; i < natoms; ++i)
        {
            const int atom_num = molinfo.number(AtomIdx(i)).value();
            auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(i));

            const int atom_idx = atom_num - 1; // 1-index versus 0-index

            Element element;

            // Infer element from mass if atomic number is negative, which
            // can be the case for files generated by Acellara's "parameterize".
            if (int(atomic_num_array[atom_idx]) < 0 and mass_array[atom_idx] > 0)
            {
                element = Element::elementWithMass(mass_array[atom_idx] * g_per_mol);
            }
            else
            {
                element = Element(int(atomic_num_array[atom_idx]));
            }

            auto amber_type = amber_type_array[atom_idx] - 1;

            params.add(AtomNum(atom_num), (charge_array[atom_idx] / AMBERCHARGECONV) * mod_electron,
                       mass_array[atom_idx] * g_per_mol, element, lj_data[amber_type],
                       ambertype_array[atom_idx].trimmed(), born_radii_array[atom_idx] * angstrom,
                       born_screening_array[atom_idx], treechains_array[atom_idx]);

            if (lj_exceptions.contains(amber_type))
                params.set(AtomNum(atom_num), lj_exceptions[amber_type]);
        }

        return params;
    };

    // function to read all of the bonding information
    auto assign_bonds = [&]()
    {
        AmberParams params(molinfo);

        const auto k_array = this->floatData("BOND_FORCE_CONSTANT").constData();
        const auto r0_array = this->floatData("BOND_EQUIL_VALUE").constData();

        auto func = [&](const QVector<int> &idxs, const QVector<qint64> &bonds, bool includes_hydrogen)
        {
            for (int idx : idxs)
            {
                const AtomNum atom0(bonds[idx] / 3 + 1);
                const AtomNum atom1(bonds[idx + 1] / 3 + 1);

                const int param_idx = bonds[idx + 2] - 1; // 1 indexed to 0 indexed

                const double k = k_array[param_idx];
                const double r0 = r0_array[param_idx];

                params.add(BondID(atom0, atom1), k, r0, includes_hydrogen);
            }
        };

        func(bonds_inc_h[molidx], this->intData("BONDS_INC_HYDROGEN"), true);
        func(bonds_exc_h[molidx], this->intData("BONDS_WITHOUT_HYDROGEN"), false);

        return params;
    };

    // function to read in the angle information
    auto assign_angles = [&]()
    {
        AmberParams params(molinfo);

        const auto k_array = this->floatData("ANGLE_FORCE_CONSTANT").constData();
        const auto t0_array = this->floatData("ANGLE_EQUIL_VALUE").constData();

        auto func = [&](const QVector<int> &idxs, const QVector<qint64> &angles, bool includes_hydrogen)
        {
            for (int idx : idxs)
            {
                const AtomNum atom0(angles[idx] / 3 + 1);
                const AtomNum atom1(angles[idx + 1] / 3 + 1);
                const AtomNum atom2(angles[idx + 2] / 3 + 1);

                const int param_idx = angles[idx + 3] - 1; // 1 indexed to 0 indexed

                const double k = k_array[param_idx];
                const double t0 = t0_array[param_idx];

                params.add(AngleID(atom0, atom1, atom2), k, t0, includes_hydrogen);
            }
        };

        func(angs_inc_h[molidx], this->intData("ANGLES_INC_HYDROGEN"), true);
        func(angs_exc_h[molidx], this->intData("ANGLES_WITHOUT_HYDROGEN"), false);

        return params;
    };

    // function to read in the dihedral+improper information and nbpairs information
    auto assign_dihedrals = [&]()
    {
        const auto k_array = this->floatData("DIHEDRAL_FORCE_CONSTANT").constData();
        const auto per_array = this->floatData("DIHEDRAL_PERIODICITY").constData();
        const auto pha_array = this->floatData("DIHEDRAL_PHASE").constData();

        const auto sceefactor = this->floatData("SCEE_SCALE_FACTOR");
        const auto scnbfactor = this->floatData("SCNB_SCALE_FACTOR");

        AmberParams params(molinfo);

        auto func = [&](const QVector<int> &idxs, const QVector<qint64> &dihedrals, bool includes_hydrogen)
        {
            for (int idx : idxs)
            {
                const AtomNum atom0(dihedrals[idx] / 3 + 1);
                const AtomNum atom1(dihedrals[idx + 1] / 3 + 1);
                const AtomNum atom2(std::abs(dihedrals[idx + 2] / 3) + 1);
                const AtomNum atom3(std::abs(dihedrals[idx + 3] / 3) + 1);

                const bool ignored = dihedrals[idx + 2] < 0;  // negative implies ignored
                const bool improper = dihedrals[idx + 3] < 0; // negative implies improper

                const int param_index = dihedrals[idx + 4] - 1; // 1 indexed to 0

                double k = k_array[param_index];     // kcal_per_mol
                double per = per_array[param_index]; // radians
                double phase = pha_array[param_index];

                if (improper)
                {
                    params.add(ImproperID(atom0, atom1, atom2, atom3), k, per, phase, includes_hydrogen);
                }
                else
                {
                    params.add(DihedralID(atom0, atom1, atom2, atom3), k, per, phase, includes_hydrogen);

                    if (not ignored)
                    {
                        // this is a 1-4 pair that has a scale factor

                        // Assume default values for 14 scaling factors
                        // Note that these are NOT inversed after reading from input
                        double sclee14 = 1.0 / AMBER14COUL;
                        double sclnb14 = 1.0 / AMBER14LJ;

                        if (not sceefactor.isEmpty())
                            sclee14 = sceefactor[param_index];

                        if (not scnbfactor.isEmpty())
                            sclnb14 = scnbfactor[param_index];

                        // invert the scale factors
                        if (sclee14 < 0.00001)
                        {
                            sclee14 = 0.0;
                        }
                        else
                        {
                            sclee14 = 1.0 / sclee14;
                        }

                        if (sclnb14 < 0.00001)
                        {
                            sclnb14 = 0.0;
                        }
                        else
                        {
                            sclnb14 = 1.0 / sclnb14;
                        }

                        // save them into the parameter space
                        params.addNB14(BondID(atom0, atom3), sclee14, sclnb14);
                    }
                }
            }
        };

        func(dihs_inc_h[molidx], this->intData("DIHEDRALS_INC_HYDROGEN"), true);
        func(dihs_exc_h[molidx], this->intData("DIHEDRALS_WITHOUT_HYDROGEN"), false);

        return params;
    };

    // function to assign all of the cmap information
    auto assign_cmaps = [&]()
    {
        AmberParams params(molinfo);

        auto func = [&](const QVector<int> &idxs, const QVector<qint64> &cmaps)
        {
            for (int idx : idxs)
            {
                // do not add 1 - these are already correctly indexed
                const AtomNum atom0(cmaps[idx]);
                const AtomNum atom1(cmaps[idx + 1]);
                const AtomNum atom2(cmaps[idx + 2]);
                const AtomNum atom3(cmaps[idx + 3]);
                const AtomNum atom4(cmaps[idx + 4]);

                const int param_idx = cmaps[idx + 5];

                params.add(atom0, atom1, atom2, atom3, atom4, cmap_data[param_idx]);
            }
        };

        func(cmap_idxs[molidx], this->intData("CMAP_INDEX") + this->intData("CHARMM_CMAP_INDEX"));

        return params;
    };

    // function to read in the excluded atoms information
    auto assign_excluded = [&]()
    {
        CLJNBPairs nbpairs;

        if (molinfo.nAtoms() <= 3)
        {
            // everything is bonded, so scale factor is 0
            nbpairs = CLJNBPairs(molinfo, CLJScaleFactor(0, 0));
        }
        else
        {
            // default is that atoms are not bonded (so scale factor is 1)
            nbpairs = CLJNBPairs(molinfo, CLJScaleFactor(1.0, 1.0));

            const int natoms = molinfo.nAtoms();

            if (molidx < 0 or molidx >= excl_atoms.count())
            {
                throw SireIO::parse_error(QObject::tr("Disagreement over the number of molecules! %1 versus %2")
                                              .arg(molidx + 1)
                                              .arg(excl_atoms.count()),
                                          CODELOC);
            }

            // get the set of excluded atoms for each atom of this molecule
            const auto excluded_atoms = excl_atoms[molidx];

            if (excluded_atoms.count() != natoms)
            {
                QStringList lines;

                for (int i = qMax(0, molidx - 3); i < qMin(molidx + 3, excl_atoms.count()); ++i)
                {
                    lines.append(QString("molidx %1: natoms = %2").arg(i).arg(excl_atoms[i].count()));
                }

                throw SireIO::parse_error(
                    QObject::tr("Disagreement over the number of atoms for molecule at index %1. "
                                "%2 versus %3. Number of atoms in neighbouring molecules are { %4 }")
                        .arg(molidx)
                        .arg(excluded_atoms.count())
                        .arg(natoms)
                        .arg(lines.join(", ")),
                    CODELOC);
            }

            for (int i = 0; i < natoms; ++i)
            {
                const AtomNum atom0 = molinfo.number(AtomIdx(i));
                const CGAtomIdx cgatom0 = molinfo.cgAtomIdx(atom0);

                // an atom is bonded to itself
                nbpairs.set(cgatom0, cgatom0, CLJScaleFactor(0, 0));

                // Are there any excluded atoms of atom0?
                for (const auto excluded_atom : excluded_atoms[i])
                {
                    if (excluded_atom > 0)
                    {
                        try
                        {
                            const AtomNum atomnum1(excluded_atom);
                            const CGAtomIdx cgatom1 = molinfo.cgAtomIdx(atomnum1);
                            nbpairs.set(cgatom0, cgatom1, CLJScaleFactor(0, 0));
                        }
                        catch (...)
                        {
                            QList<int> nums;

                            for (int iat = 0; iat < molinfo.nAtoms(); ++iat)
                            {
                                nums.append(molinfo.number(AtomIdx(iat)));
                            }

                            throw SireIO::parse_error(
                                QObject::tr("Cannot find atom %1 in the molecule. This is likely because "
                                            "the code to deal with molecules with SS bonds is broken..."
                                            "Valid atom numbers are %2. Listed excluded atoms are %3")
                                    .arg(excluded_atom)
                                    .arg(Sire::toString(nums))
                                    .arg(Sire::toString(excluded_atoms[i])),
                                CODELOC);
                        }
                    }
                }
            }
        }

        AmberParams params(molinfo);
        params.setExcludedAtoms(nbpairs);

        return params;
    };

    AmberParams params(molinfo);

    // assign all of the parameters
    if (usesParallel())
    {
        AmberParams atoms, bonds, angles, dihedrals, excluded, cmaps;

        tbb::parallel_invoke([&]()
                             { atoms = assign_atoms(); },
                             [&]()
                             { bonds = assign_bonds(); },
                             [&]()
                             { angles = assign_angles(); },
                             [&]()
                             { dihedrals = assign_dihedrals(); },
                             [&]()
                             { excluded = assign_excluded(); },
                             [&]()
                             { cmaps = assign_cmaps(); });

        params += atoms;
        params += bonds;
        params += angles;
        params += dihedrals;
        params += excluded;
        params += cmaps;
    }
    else
    {
        params += assign_atoms();
        params += assign_bonds();
        params += assign_angles();
        params += assign_dihedrals();
        params += assign_excluded();
        params += assign_cmaps();
    }

    const auto radius_set = this->stringData("RADIUS_SET");

    if (radius_set.count() > 0)
    {
        params.setRadiusSet(radius_set[0]);
    }

    return params;
}

/** Internal function to set a property in a molecule */
void _setProperty(MolEditor &mol, const PropertyMap &map, QString key, const Property &value)
{
    const auto mapped = map[key];

    if (mapped.hasValue())
    {
        mol.setProperty(key, mapped.value());
    }
    else
    {
        mol.setProperty(mapped.source(), value);
    }
}

/** Internal function used to get the molecule structure that starts at index 'start_idx'
    in the file, and that has 'natoms' atoms */
MolEditor AmberPrm::getMoleculeEditor(int molidx, const PropertyMap &map) const
{
    // first, construct the layout of the molecule (sorting of atoms into residues and cutgroups)
    auto mol = this->getMolStructure(molidx, map["cutting"]);

    if (mol.nAtoms() == 0)
        return MolEditor();

    // get the info object that can map between AtomNum to AtomIdx etc.
    const auto molinfo = mol.info();

    auto amber_params = this->getAmberParams(molidx, molinfo);

    amber_params.setPropertyMap(map);

    _setProperty(mol, map, "charge", amber_params.charges());
    _setProperty(mol, map, "LJ", amber_params.ljs());
    _setProperty(mol, map, "mass", amber_params.masses());
    _setProperty(mol, map, "element", amber_params.elements());
    _setProperty(mol, map, "ambertype", amber_params.amberTypes());
    _setProperty(mol, map, "atomtype", amber_params.amberTypes());
    _setProperty(mol, map, "connectivity", amber_params.connectivity());
    _setProperty(mol, map, "bond", amber_params.bondFunctions(InternalPotential::symbols().bond().r()));
    _setProperty(mol, map, "angle", amber_params.angleFunctions(InternalPotential::symbols().angle().theta()));
    _setProperty(mol, map, "dihedral", amber_params.dihedralFunctions(InternalPotential::symbols().dihedral().phi()));
    _setProperty(mol, map, "improper", amber_params.improperFunctions(InternalPotential::symbols().dihedral().phi()));
    _setProperty(mol, map, "intrascale", amber_params.cljScaleFactors());
    _setProperty(mol, map, "gb_radii", amber_params.gbRadii());
    _setProperty(mol, map, "gb_screening", amber_params.gbScreening());
    _setProperty(mol, map, "gb_radius_set", StringProperty(amber_params.radiusSet()));
    _setProperty(mol, map, "treechain", amber_params.treeChains());
    _setProperty(mol, map, "parameters", amber_params);
    _setProperty(mol, map, "forcefield", ffield);

    if (not amber_params.cmapFunctions().isEmpty())
    {
        _setProperty(mol, map, "cmap", amber_params.cmapFunctions());
    }

    return mol;
}

/** Return the ith molecule that is described by this AmberPrm file. Note
    that this molecule won't have any coordinate data, as this is not
    provided in this file */
Molecule AmberPrm::getMolecule(int molidx, const PropertyMap &map) const
{
    molidx = Index(molidx).map(this->nMolecules());
    return this->getMoleculeEditor(molidx, map).commit();
}

/** Return the ith molecule that is described by this AmberPrm file, getting
    the coordinate (and possibly velocity) data from the passed AmberRst file */
Molecule AmberPrm::getMolecule(int molidx, const AmberRst7 &rst, const PropertyMap &map) const
{
    if (rst.nAtoms() != this->nAtoms())
    {
        // these two files are likely to be incompatible!
        throw SireIO::parse_error(QObject::tr("The passed Amber Parm and Restart files are incompatible as they "
                                              "contain data for different numbers of atoms "
                                              "(%1 in the Parm file and %2 in the Restart)")
                                      .arg(this->nAtoms())
                                      .arg(rst.nAtoms()),
                                  CODELOC);
    }

    molidx = Index(molidx).map(this->nMolecules());

    auto mol = this->getMoleculeEditor(molidx, map);

    if (mol.nAtoms() == 0)
        return mol;

    const auto molinfo = mol.info();

    // create space for the coordinates
    auto coords = QVector<QVector<Vector>>(molinfo.nCutGroups());

    for (int i = 0; i < molinfo.nCutGroups(); ++i)
    {
        coords[i] = QVector<Vector>(molinfo.nAtoms(CGIdx(i)));
    }

    const Vector *coords_array = rst.coordinates().constData();

    if (rst.hasVelocities())
    {
        auto vels = AtomVelocities(molinfo);
        const Vector *vels_array = rst.velocities().constData();

        for (int i = 0; i < mol.nAtoms(); ++i)
        {
            const int atom_num = molinfo.number(AtomIdx(i)).value();
            auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(i));

            const int atom_idx = atom_num - 1; // 1-index versus 0-index

            coords[cgatomidx.cutGroup()][cgatomidx.atom()] = coords_array[atom_idx];

            // velocity is Angstroms per 1/20.455 ps
            const auto vel_unit = (1.0 / 20.455) * angstrom / picosecond;

            const Vector &vel = vels_array[atom_idx];
            vels.set(cgatomidx, Velocity3D(vel.x() * vel_unit, vel.y() * vel_unit, vel.z() * vel_unit));
        }

        return mol.setProperty(map["velocity"], vels)
            .setProperty(map["coordinates"], AtomCoords(CoordGroupArray(coords)))
            .commit();
    }
    else
    {
        for (int i = 0; i < mol.nAtoms(); ++i)
        {
            const int atom_num = molinfo.number(AtomIdx(i)).value();
            auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(i));

            const int atom_idx = atom_num - 1; // 1-index versus 0-index

            coords[cgatomidx.cutGroup()][cgatomidx.atom()] = coords_array[atom_idx];
        }

        return mol.setProperty(map["coordinates"], AtomCoords(CoordGroupArray(coords))).commit();
    }
}

/** Return the format name that is used to identify this file format within Sire */
QString AmberPrm::formatName() const
{
    return "PRM7";
}

/** Return the suffixes that AmberPrm files are normally associated with */
QStringList AmberPrm::formatSuffix() const
{
    static const QStringList suffixes = {"prm7", "prm", "top7", "top", "prmtop7", "prmtop"};
    return suffixes;
}

bool AmberPrm::isTopology() const
{
    return true;
}

/** Return a description of the file format */
QString AmberPrm::formatDescription() const
{
    return QObject::tr("Amber topology/parameter format files supported from Amber 7 upwards.");
}

/** Return the System that is described by this AmberPrm file. Note that
    the molecules in this system don't have any coordinates (as these aren't
    provided by the file */
System AmberPrm::startSystem(const PropertyMap &map) const
{
    const int nmols = this->nMolecules();

    if (nmols == 0)
        return System();

    QVector<MolEditor> mols(nmols);
    auto *mols_array = mols.data();

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                          {
            // create and populate all of the molecules
            for (int i = r.begin(); i < r.end(); ++i)
            {
                mols_array[i] = this->getMoleculeEditor(i, map);
            } });
    }
    else
    {
        for (int i = 0; i < nmols; ++i)
        {
            mols_array[i] = this->getMoleculeEditor(i, map);
        }
    }

    // we need to check that the atoms are still in the same order
    // as they appeared in the prmtop file - the order can change if
    // molecules are discontiguous (i.e. ions moved to the end but
    // bonded to earlier chains of a protein)
    int expected_atomnum = 1;
    bool in_expected_order = true;
    QVector<qint64> loaded_order;

    ProgressBar bar("Checking atom order");
    bool entered_bar = false;

    for (auto &mol : mols)
    {
        for (int i = 0; i < mol.nAtoms(); ++i)
        {
            const auto atomnum = mol.info().number(AtomIdx(i));

            loaded_order.append(atomnum.value() - 1);

            if (atomnum.value() != expected_atomnum)
            {
                if (not entered_bar)
                {
                    entered_bar = true;
                    bar.setSpeedUnit("atoms / s");
                    bar = bar.enter();
                }

                in_expected_order = false;

                // we need to renumber this atom so that it has the
                // expected atom number - they must increase in order
                // according to molidx/atomidx order
                bar.tick(QString("Renumbering atom %1 to %2").arg(atomnum.value()).arg(expected_atomnum));

                auto atom = mol.atom(AtomIdx(i));
                atom.renumber(AtomNum(expected_atomnum));
                mol = atom.molecule();
            }

            expected_atomnum += 1;
        }
    }

    if (entered_bar)
    {
        bar.success();
    }

    MoleculeGroup molgroup("all");

    for (auto mol : mols)
    {
        if (mol.nAtoms() > 0)
            molgroup.add(mol.commit());
    }

    System system(this->title());
    system.add(molgroup);
    system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));

    if (not in_expected_order)
    {
        system.setProperty("loaded_atom_order", IntegerArrayProperty(loaded_order));
    }

    // some top files contains "BOX_DIMENSIONS" information. Add this now, as it
    // now in case it is not replaced by the coordinates file
    const auto boxdims = float_data.value("BOX_DIMENSIONS");

    if (boxdims.count() == 4)
    {
        // ANGLE, X, Y, Z dimensions - we will only read this if the angle is 90
        if (boxdims[0] == 90.0)
        {
            PeriodicBox space(Vector(boxdims[1], boxdims[2], boxdims[3]));
            system.setProperty("space", space);
        }
    }

    return system;
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr AmberPrm::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberPrm(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr AmberPrm::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberPrm(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr AmberPrm::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberPrm(system, map));
}
