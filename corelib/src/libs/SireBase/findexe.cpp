/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "findexe.h"

#include <QDir>
#include <QProcess>
#include <QRegExp>
#include <QStringList>

#include "SireError/errors.h"

QFileInfo SireBase::findExe(const QString &exe)
{
#ifdef Q_OS_WIN
    const QString splitChar(";");
    const QString exeSuffix = ".exe";
#else
    const QString splitChar(":");
    const QString exeSuffix;
#endif
    // does it exist in the current directory?
    QFileInfo progfile(exe);

    if (progfile.isExecutable())
        return progfile;

    // try to find exe in the application PATH
    QStringList env_variables = QProcess::systemEnvironment();

    QRegExp rexp("^PATH=(.*)");
    rexp.setCaseSensitivity(Qt::CaseInsensitive);

    QStringList path;

    foreach (QString env_variable, env_variables)
    {
        if (rexp.indexIn(env_variable) != -1)
        {
            path = rexp.cap(1).split(splitChar, Qt::SkipEmptyParts);
            break;
        }
    }

    foreach (QString dir, path)
    {
        QFileInfo new_exe(QDir(dir).filePath(progfile.fileName() + exeSuffix));

        if (new_exe.isExecutable())
        {
            // we have found the executable!
            return new_exe;
        }
    }

    throw SireError::process_error(QObject::tr("Could not find the executable \"%1\" in the current directory, nor "
                                               "in the system PATH.")
                                       .arg(exe),
                                   CODELOC);

    return QFileInfo();
}
