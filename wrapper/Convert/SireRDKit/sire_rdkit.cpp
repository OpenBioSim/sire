
#include "sire_rdkit.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MonomerInfo.h>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atommasses.h"
#include "SireMol/atommatch.h"
#include "SireMol/atomproperty.hpp"
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"
#include "SireMol/bondorder.h"
#include "SireMol/stereochemistry.h"
#include "SireMol/chirality.h"
#include "SireMol/hybridization.h"
#include "SireMol/iswater.h"
#include "SireMol/mover_metaid.h"

#include "SireMM/selectorbond.h"

#include "SireMaths/vector.h"

#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"
#include "SireBase/stringproperty.h"

#include "SireUnits/units.h"

#include "tostring.h"

#include <map>

#include <QDebug>

using RDKit::ROMOL_SPTR;
using RDKit::RWMOL_SPTR;
using SireBase::PropertyMap;
using SireMol::Molecule;
using SireMol::SelectorMol;

namespace SireRDKit
{
    bool use_parallel(int n, const SireBase::PropertyMap &map)
    {
#ifdef Q_OS_WIN
        // parallel RDKit segfaults on windows!
        return false;
#else
        if (n <= 16)
            return false;

        if (map["parallel"].hasValue())
        {
            return map["parallel"].value().asABoolean();
        }

        return true;
#endif
    }

    RDKit::Atom::ChiralType string_to_chiral(const QString &typ)
    {
        const static std::map<QString, RDKit::Atom::ChiralType> types = {
            {"CHI_UNSPECIFIED", RDKit::Atom::CHI_UNSPECIFIED},
            {"CHI_TETRAHEDRAL_CW", RDKit::Atom::CHI_TETRAHEDRAL_CW},
            {"CHI_TETRAHEDRAL_CCW", RDKit::Atom::CHI_TETRAHEDRAL_CCW},
            {"CHI_OTHER", RDKit::Atom::CHI_OTHER},
            // These don't seem to be supported on older RDKits?
            //  {"CHI_TETRAHEDRAL", RDKit::Atom::CHI_OTHER},
            //  {"CHI_ALLENE", RDKit::Atom::CHI_ALLENE},
            //  {"CHI_SQUAREPLANAR", RDKit::Atom::CHI_SQUAREPLANAR},
            //  {"CHI_TRIGONALBIPYRAMIDAL", RDKit::Atom::CHI_TRIGONALBIPYRAMIDAL},
            //  {"CHI_OCTAHEDRAL", RDKit::Atom::CHI_OCTAHEDRAL},
            {"CHI_UNKNOWN", RDKit::Atom::CHI_UNSPECIFIED}};

        auto it = types.find(typ);

        if (it == types.end())
            return RDKit::Atom::CHI_UNSPECIFIED;
        else
            return it->second;
    }

    QString chiral_to_string(RDKit::Atom::ChiralType typ)
    {
        switch (typ)
        {
        case RDKit::Atom::CHI_UNSPECIFIED:
            return "CHI_UNSPECIFIED";
        case RDKit::Atom::CHI_TETRAHEDRAL_CW:
            return "CHI_TETRAHEDRAL_CW";
        case RDKit::Atom::CHI_TETRAHEDRAL_CCW:
            return "CHI_TETRAHEDRAL_CCW";
        case RDKit::Atom::CHI_OTHER:
            return "CHI_OTHER";
            // These don't seem to be supported on older RDKits?
            //   case RDKit::Atom::CHI_TETRAHEDRAL:
            //       return "CHI_TETRAHEDRAL";
            //   case RDKit::Atom::CHI_ALLENE:
            //       return "CHI_ALLENE";
            //   case RDKit::Atom::CHI_SQUAREPLANAR:
            //       return "CHI_SQUAREPLANAR";
            //   case RDKit::Atom::CHI_TRIGONALBIPYRAMIDAL:
            //       return "CHI_TRIGONALBIPYRAMIDAL";
            //   case RDKit::Atom::CHI_OCTAHEDRAL:
            //       return "CHI_OCTAHEDRAL";
        default:
            return "CHI_UNKNOWN";
        }
    }

    RDKit::Atom::HybridizationType string_to_hybridization(const QString &typ)
    {
        const static std::map<QString, RDKit::Atom::HybridizationType> types = {
            {"UNSPECIFIED", RDKit::Atom::UNSPECIFIED},
            {"S", RDKit::Atom::S},
            {"SP", RDKit::Atom::SP},
            {"SP2", RDKit::Atom::SP2},
            {"SP3", RDKit::Atom::SP3},
            // These don't seem to be supported on older RDKits?
            // {"SP2D", RDKit::Atom::SP2D},
            {"SP3D", RDKit::Atom::SP3D},
            {"SP3D2", RDKit::Atom::SP3D2},
            {"OTHER", RDKit::Atom::OTHER},
            {"UNKNOWN", RDKit::Atom::UNSPECIFIED}};

        auto it = types.find(typ);

        if (it == types.end())
            return RDKit::Atom::UNSPECIFIED;
        else
            return it->second;
    }

    QString hybridization_to_string(RDKit::Atom::HybridizationType typ)
    {
        switch (typ)
        {
        case RDKit::Atom::UNSPECIFIED:
            return "UNSPECIFIED";
        case RDKit::Atom::S:
            return "S";
        case RDKit::Atom::SP:
            return "SP";
        case RDKit::Atom::SP2:
            return "SP2";
        case RDKit::Atom::SP3:
            return "SP3";
            // These don't seem to be supported on older RDKits?
            //    case RDKit::Atom::SP2D:
            //        return "SP2D";
        case RDKit::Atom::SP3D:
            return "SP3D";
        case RDKit::Atom::SP3D2:
            return "SP3D2";
        case RDKit::Atom::OTHER:
            return "OTHER";
        default:
            return "UNKNOWN";
        }
    }

    QString stereo_to_string(RDKit::Bond::BondStereo stereo)
    {
        switch (stereo)
        {
        case RDKit::Bond::STEREONONE:
            return "STEREONONE";
        case RDKit::Bond::STEREOANY:
            return "STEREOANY";
        case RDKit::Bond::STEREOZ:
            return "STEREOZ";
        case RDKit::Bond::STEREOE:
            return "STEREOE";
        case RDKit::Bond::STEREOCIS:
            return "STEREOCIS";
        case RDKit::Bond::STEREOTRANS:
            return "STEREOTRANS";
        default:
            return "STEREONONE";
        }
    }

    RDKit::Bond::BondStereo string_to_stereo(const QString &stereo)
    {
        const static std::map<QString, RDKit::Bond::BondStereo> stereos = {
            {"STEREONONE", RDKit::Bond::STEREONONE},
            {"STEREOANY", RDKit::Bond::STEREOANY},
            {"STEREOZ", RDKit::Bond::STEREOZ},
            {"STEREOE", RDKit::Bond::STEREOE},
            {"STEREOCIS", RDKit::Bond::STEREOCIS},
            {"STEREOTRANS", RDKit::Bond::STEREOTRANS}};

        auto it = stereos.find(stereo);

        if (it == stereos.end())
            return RDKit::Bond::STEREONONE;
        else
            return it->second;
    }

    QString bondtype_to_string(RDKit::Bond::BondType typ)
    {
        switch (typ)
        {
        case RDKit::Bond::UNSPECIFIED:
            return "UNSPECIFIED";
        case RDKit::Bond::SINGLE:
            return "SINGLE";
        case RDKit::Bond::DOUBLE:
            return "DOUBLE";
        case RDKit::Bond::TRIPLE:
            return "TRIPLE";
        case RDKit::Bond::QUADRUPLE:
            return "QUADRUPLE";
        case RDKit::Bond::QUINTUPLE:
            return "QUINTUPLE";
        case RDKit::Bond::HEXTUPLE:
            return "HEXTUPLE";
        case RDKit::Bond::ONEANDAHALF:
            return "ONEANDAHALF";
        case RDKit::Bond::TWOANDAHALF:
            return "TWOANDAHALF";
        case RDKit::Bond::THREEANDAHALF:
            return "THREEANDAHALF";
        case RDKit::Bond::FOURANDAHALF:
            return "FOURANDAHALF";
        case RDKit::Bond::FIVEANDAHALF:
            return "FIVEANDAHALF";
        case RDKit::Bond::AROMATIC:
            return "AROMATIC";
        case RDKit::Bond::IONIC:
            return "IONIC";
        case RDKit::Bond::HYDROGEN:
            return "HYDROGEN";
        case RDKit::Bond::THREECENTER:
            return "THREECENTER";
        case RDKit::Bond::DATIVEONE:
            return "DATIVEONE";
        case RDKit::Bond::DATIVE:
            return "DATIVEONE";
        case RDKit::Bond::DATIVEL:
            return "DATIVEL";
        case RDKit::Bond::DATIVER:
            return "DATIVER";
        case RDKit::Bond::OTHER:
            return "OTHER";
        case RDKit::Bond::ZERO:
            return "ZERO";
        default:
            return "UNSPECIFIED";
        }
    }

    RDKit::Bond::BondType order_to_bondtype(int order)
    {
        switch (order)
        {
        case 0:
            return RDKit::Bond::ZERO;
        case 1:
            return RDKit::Bond::SINGLE;
        case 2:
            return RDKit::Bond::DOUBLE;
        case 3:
            return RDKit::Bond::TRIPLE;
        case 4:
            return RDKit::Bond::QUADRUPLE;
        case 5:
            return RDKit::Bond::HEXTUPLE;
        default:
            return RDKit::Bond::OTHER;
        }
    }

    RDKit::Bond::BondType string_to_bondtype(const QString &typ)
    {
        const static std::map<QString, RDKit::Bond::BondType> types = {
            {"UNSPECIFIED", RDKit::Bond::UNSPECIFIED},
            {"SINGLE", RDKit::Bond::SINGLE},
            {"DOUBLE", RDKit::Bond::DOUBLE},
            {"TRIPLE", RDKit::Bond::TRIPLE},
            {"QUADRUPLE", RDKit::Bond::QUADRUPLE},
            {"QUINTUPLE", RDKit::Bond::QUINTUPLE},
            {"HEXTUPLE", RDKit::Bond::HEXTUPLE},
            {"ONEANDAHALF", RDKit::Bond::ONEANDAHALF},
            {"TWOANDAHALF", RDKit::Bond::TWOANDAHALF},
            {"THREEANDAHALF", RDKit::Bond::THREEANDAHALF},
            {"FOURANDAHALF", RDKit::Bond::FOURANDAHALF},
            {"FIVEANDAHALF", RDKit::Bond::FIVEANDAHALF},
            {"AROMATIC", RDKit::Bond::AROMATIC},
            {"IONIC", RDKit::Bond::IONIC},
            {"HYDROGEN", RDKit::Bond::HYDROGEN},
            {"THREECENTER", RDKit::Bond::THREECENTER},
            {"DATIVEONE", RDKit::Bond::DATIVEONE},
            {"DATIVE", RDKit::Bond::DATIVE},
            {"DATIVEL", RDKit::Bond::DATIVEL},
            {"DATIVER", RDKit::Bond::DATIVER},
            {"OTHER", RDKit::Bond::OTHER},
            {"ZERO", RDKit::Bond::ZERO}};

        auto it = types.find(typ);

        if (it == types.end())
            return RDKit::Bond::UNSPECIFIED;
        else
            return it->second;
    }

    template <class T, class V>
    void set_prop(T &atom, const QString &key, const V &value,
                  const PropertyMap &map)
    {
        const auto name = map[key];

        if (name.hasSource())
            atom.setProperty(name.source(), value);
    }

    /**
     *
        """Calculate the number of unpaired electrons (NUE) of an atom

        Parameters
        ----------
        atom: rdkit.Chem.rdchem.Atom
            The atom for which the NUE will be computed

        Returns
        -------
        nue : list
            The NUE for each possible valence of the atom
        """

    */
    QVector<int> get_nb_unpaired_electrons(const RDKit::Atom &atom)
    {
        const RDKit::PeriodicTable *PERIODIC_TABLE = RDKit::PeriodicTable::getTable();

        auto expected_vs = PERIODIC_TABLE->getValenceList(atom.getAtomicNum());

        auto current_v = atom.getTotalValence() - atom.getFormalCharge();

        QVector<int> ret;
        ret.reserve(expected_vs.size());

        for (int i = 0; i < expected_vs.size(); ++i)
        {
            ret.append(expected_vs[i] - current_v);
        }

        return ret;
    }

    /** Infer the bond info and charge state for all atoms of 'molecule'.
     *  This is heavily based on _infer_bond_info from MDAnalysis,
     *  as described here:
     *  https://blog.matteoferla.com/2020/02/guess-bond-order-in-rdkit-by-number-of.html
     *
     *  https://docs.mdanalysis.org/2.0.0/_modules/MDAnalysis/converters/RDKit.html#_infer_bo_and_charges
     *
     *  MDAnalysis code was written by Matteo Ferla and was released
     *  under the GPL. This is a C++ port of the original python code
     *
     *  The molecule is edited in place and must have all hydrogens added.
     *  Quoting the MDAnalysis docs
     *
        Since most MD topology files don't explicitly retain information on bond
        orders or charges, it has to be guessed from the topology. This is done by
        looping over each atom and comparing its expected valence to the current
        valence to get the Number of Unpaired Electrons (NUE).
        If an atom has a negative NUE, it needs a positive formal charge (-NUE).
        If two neighbouring atoms have UEs, the bond between them most
        likely has to be increased by the value of the smallest NUE.
        If after this process, an atom still has UEs, it needs a negative formal
        charge of -NUE.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.RWMol
            The molecule is modified inplace and must have all hydrogens added

        Notes
        -----
        This algorithm is order dependant. For example, for a carboxylate group
        R-C(-O)-O the first oxygen read will receive a double bond and the other
        one will be charged. It will also affect more complex conjugated systems.
     *
     */
    void infer_bond_info(RDKit::RWMol &molecule)
    {
        // get all the atoms
        QList<std::pair<QVector<int>, RDKit::Atom *>> atoms;

        for (auto atom : molecule.atoms())
        {
            if (atom->getAtomicNum() > 1)
            {
                atom->setNoImplicit(true);
                atoms.append(std::make_pair(get_nb_unpaired_electrons(*atom),
                                            atom));
            }
        }

        // sort these atoms so that the ones with most unpaired electrons
        // come first
        std::stable_sort(atoms.begin(), atoms.end(),
                         [](const std::pair<QVector<int>, RDKit::Atom *> &atom0,
                            const std::pair<QVector<int>, RDKit::Atom *> &atom1)
                         {
                             return std::get<0>(atom0).at(0) > std::get<0>(atom1).at(0);
                         });

        for (const auto &p : atoms)
        {
            // for this atom...
            RDKit::Atom *atom = std::get<1>(p);

            if (atom->getDegree() == 0)
            {
                // no neighbors, so no bonds - this could be a monovalent cation
                switch (atom->getAtomicNum())
                {
                case 3:
                case 11:
                case 19:
                case 37:
                case 47:
                case 55:
                    atom->setFormalCharge(1);
                    break;
                case 12:
                case 20:
                case 29:
                case 30:
                case 38:
                case 56:
                case 26:
                    atom->setFormalCharge(2);
                    break;
                case 13:
                    // Fe could also be + 3
                    atom->setFormalCharge(3);
                    break;
                default:
                    // no, it is an anion - use the negative of the number of
                    // unpaired electrons
                    atom->setFormalCharge(-get_nb_unpaired_electrons(*atom)[0]);
                    break;
                }

                molecule.updatePropertyCache(false);
                continue;
            }

            // number of unpaired electrons
            auto nue = get_nb_unpaired_electrons(*atom);

            // if there's only one possible valence state and the corresponding
            // nue is negative, it means we can only add a positive charge to
            // the atom
            if (nue.count() == 1 and nue[0] < 0)
            {
                atom->setFormalCharge(-nue[0]);
                molecule.updatePropertyCache(false);
            }

            if (nue.count() == 1 and nue[0] <= 0)
            {
                // go to the next atom if above case or atom has
                // no unpaired electrons
                continue;
            }

            // get the neighbors
            QList<std::pair<QVector<int>, RDKit::Atom *>> neighbors;

            for (auto n : molecule.atomNeighbors(atom))
            {
                neighbors.append(std::make_pair(get_nb_unpaired_electrons(*n),
                                                n));
            }

            std::stable_sort(neighbors.begin(), neighbors.end(),
                             [](const std::pair<QVector<int>, RDKit::Atom *> &atom0,
                                const std::pair<QVector<int>, RDKit::Atom *> &atom1)
                             {
                                 return std::get<0>(atom0).at(0) > std::get<0>(atom1).at(0);
                             });

            for (int i = 0; i < neighbors.count(); ++i)
            {
                int min_nue = 0;

                for (const auto n : nue)
                {
                    if (n > 0)
                    {
                        if (min_nue == 0 or n < min_nue)
                            min_nue = n;
                    }
                }

                const auto na_nue = std::get<0>(neighbors.at(i));
                auto neighbor = std::get<1>(neighbors.at(i));

                int min_na_nue = 0;

                for (const auto n : na_nue)
                {
                    if (n > 0)
                    {
                        if (min_na_nue == 0 or n < min_na_nue)
                            min_na_nue = n;
                    }
                }

                // check if one of the neighbors has a common nue
                auto common_nue = std::min(min_nue, min_na_nue);

                // a common nue of 0 means we don't need to do anything
                if (common_nue != 0)
                {
                    auto bond = molecule.getBondBetweenAtoms(
                        atom->getIdx(), neighbor->getIdx());

                    int order = common_nue + 1;
                    bond->setBondType(order_to_bondtype(order));
                    molecule.updatePropertyCache(false);

                    if (i < neighbors.count())
                    {
                        // recalculate the nue for this atom
                        nue = get_nb_unpaired_electrons(*atom);
                    }

                    bool any_equal_zero = false;

                    for (auto n : nue)
                    {
                        if (n == 0)
                        {
                            any_equal_zero = true;
                            break;
                        }
                    }

                    if (any_equal_zero)
                    {
                        // we can stop here
                        break;
                    }
                }
            } // end of loop over neighbours

            // recalculate the nue for this atom
            nue = get_nb_unpaired_electrons(*atom);

            bool any_equal_zero = false;

            for (auto n : nue)
            {
                if (n == 0)
                {
                    any_equal_zero = true;
                    break;
                }
            }

            if (not any_equal_zero)
            {
                // transform it to a negative charge
                atom->setFormalCharge(-nue[0]);
                atom->setNumRadicalElectrons(0);
                molecule.updatePropertyCache(false);
            }
        }
    }

    ROMOL_SPTR sire_to_rdkit(const Molecule &mol, const PropertyMap &map)
    {
        RDKit::RWMol molecule;
        molecule.beginBatchEdit();

        // set the name of the molecule
        std::string name;
        if (mol.hasProperty("name"))
        {
            name = mol.property("name").asAString().toStdString();
        }
        else
        {
            name = mol.name().value().toStdString();
        }
        molecule.setProp<std::string>("_Name", name);

        // set any SDF tags as properties
        std::string sdf_tag;
        if (mol.hasProperty("sdf_data"))
        {
            const auto sdf_data = mol.property("sdf_data").asA<SireBase::Properties>();

            for (const auto &tag : sdf_data.propertyKeys())
            {
                try
                {
                    molecule.setProp<std::string>(tag.toStdString(), sdf_data.property(tag).asAString().toStdString());
                }
                catch (...)
                {
                    const auto string_array = sdf_data.property(tag).asA<SireBase::StringArrayProperty>();

                    QString string;
                    for (int i=0; i<string_array.size(); i++)
                    {
                        string.append(string_array[i] + "\n");
                    }

                    molecule.setProp<std::string>(tag.toStdString(), string.toStdString());
                }
            }
        }

        const auto atoms = mol.atoms();

        QList<SireMol::Element> elements;

        QVector<SireMaths::Vector> coords(atoms.count());
        SireMaths::Vector *coords_data = coords.data();
        bool has_coords = false;

        bool force_stereo_inference = false;
        if (map.specified("force_stereo_inference"))
        {
            force_stereo_inference = map["force_stereo_inference"].value().asABoolean();
        }

        for (int i = 0; i < atoms.count(); ++i)
        {
            const auto atom = atoms(i);

            molecule.addAtom(true);
            auto a = molecule.getActiveAtom();

            // create a AtomPDBResidueInfo object for the atom, and
            // populate it with the name and residue information
            auto info = new RDKit::AtomPDBResidueInfo();

            info->setSerialNumber(atom.number().value());
            info->setName(atom.name().value().toStdString());

            if (atom.isWithinResidue())
            {
                auto residue = atom.residue();
                info->setResidueName(residue.name().value().toStdString());
                info->setResidueNumber(residue.number().value());
            }

            if (atom.isWithinChain())
            {
                auto chain = atom.chain();
                info->setChainId(chain.name().value().toStdString());
            }

            a->setMonomerInfo(info);
            a->setProp("molFileAlias", info->getName());

            const auto element = atom.property<SireMol::Element>(map["element"]);

            a->setAtomicNum(element.nProtons());

            elements.append(element);

            try
            {
                coords_data[i] = atom.property<SireMaths::Vector>(map["coordinates"]);
                has_coords = true;
            }
            catch (...)
            {
            }

            try
            {
                a->setFormalCharge(atom.property<SireUnits::Dimension::Charge>(map["formal_charge"]).to(SireUnits::mod_electron));
            }
            catch (...)
            {
            }

            try
            {
                a->setIsAromatic(atom.property<qint64>(map["is_aromatic"]) != 0);
            }
            catch (...)
            {
            }

            try
            {
                a->setIsotope(atom.property<qint64>(map["isotope"]));
            }
            catch (...)
            {
            }

            try
            {
                a->setChiralTag(string_to_chiral(atom.property<SireMol::Chirality>(map["chirality"]).toRDKit()));
            }
            catch (...)
            {
            }

            try
            {
                a->setHybridization(string_to_hybridization(atom.property<SireMol::Hybridization>(map["hybridization"]).toRDKit()));
            }
            catch (...)
            {
            }
        }

        const auto bonds = SireMM::SelectorBond(mol, map);

        bool has_bond_info = false;

        for (int i = 0; i < bonds.count(); ++i)
        {
            const auto bond = bonds(i);

            RDKit::Bond::BondType bondtype = RDKit::Bond::SINGLE;
            RDKit::Bond::BondStereo stereo = RDKit::Bond::STEREONONE;

            try
            {
                bondtype = string_to_bondtype(bond.property(map["order"]).asA<SireMol::BondOrder>().toRDKit());

                // one bond has bond info, so assume that all do
                has_bond_info = true;
            }
            catch (...)
            {
            }

            try
            {
                stereo = string_to_stereo(bond.property(map["stereochemistry"]).asA<SireMol::Stereochemistry>().toRDKit());
            }
            catch (...)
            {
            }

            const auto atom0 = bond.atom0().index().value();
            const auto atom1 = bond.atom1().index().value();

            if (elements.at(atom0).nProtons() == 0 or
                elements.at(atom1).nProtons() == 0)
            {
                // bonds involving dummy atoms
                bondtype = RDKit::Bond::ZERO;
                stereo = RDKit::Bond::STEREONONE;
            }
            else if (elements.at(atom0).nProtons() == 1 and
                     elements.at(atom1).nProtons() == 1 and
                     bonds.count() > 1)
            {
                // H-H bonds in molecules containing more
                // bonds than just hydrogens...
                bondtype = RDKit::Bond::ZERO;
                stereo = RDKit::Bond::STEREONONE;
            }

            molecule.addBond(bond.atom0().index().value(),
                             bond.atom1().index().value(),
                             bondtype);
            molecule.getBondWithIdx(i)->setStereo(stereo);
        }

        if (has_coords)
        {
            const int nats = coords.count();
            const SireMaths::Vector *data = coords.constData();

            RDKit::Conformer *conformer = new RDKit::Conformer(nats);

            try
            {
                for (int i = 0; i < nats; ++i)
                {
                    const SireMaths::Vector &p = data[i];
                    conformer->setAtomPos(i, RDGeom::Point3D(p.x(), p.y(), p.z()));
                }

                // the molecule should take over ownership of this pointer now
                molecule.addConformer(conformer, true);
            }
            catch (...)
            {
                delete conformer;
                throw;
            }
        }

        molecule.commitBatchEdit();
        molecule.updatePropertyCache(false);

        if (atoms.count() > 1 and (not has_bond_info or force_stereo_inference))
        {
            // we need to infer the bond information
            infer_bond_info(molecule);
        }

        // sanitze the molecule.
        try
        {
            RDKit::MolOps::sanitizeMol(molecule);
        }
        catch (...)
        {
        }

        // try assigning stereochemistry from 3D coordinates as it is the most
        // reliable way to do it
        if (has_coords and not force_stereo_inference)
        {
            try
            {
                RDKit::MolOps::assignStereochemistryFrom3D(molecule);
            }
            catch (...)
            {
            }
        }

        return ROMOL_SPTR(new RDKit::ROMol(molecule));
    }

    Molecule rdkit_to_sire(const ROMOL_SPTR &mol, const PropertyMap &map)
    {
        if (mol.get() == 0)
            return Molecule();

        QString molname;

        if (mol->hasProp("_Name"))
        {
            molname = QString::fromStdString(mol->getProp<std::string>("_Name"));
        }

        auto cg = Molecule().edit().rename(molname).add(SireMol::CGName("0"));

        int n = 0;

        QHash<const RDKit::Atom *, SireMol::AtomIdx> atom_to_idx;

        QHash<SireMol::Element, int> element_count;

        for (const auto &atom : mol->atoms())
        {
            atom_to_idx[atom] = n;
            n += 1;

            // see if there is a PDBResidueInfo object for this atom
            const auto *info = atom->getMonomerInfo();

            SireMol::AtomName atomname;
            SireMol::AtomNum atomnum;
            SireMol::ResName resname;
            SireMol::ResNum resnum;
            SireMol::ChainName chainname;

            if (info != 0)
            {
                atomname = SireMol::AtomName(QString::fromStdString(info->getName()));

                if (info->getMonomerType() == RDKit::AtomMonomerInfo::PDBRESIDUE)
                {
                    auto resinfo = static_cast<const RDKit::AtomPDBResidueInfo *>(info);

                    atomname = SireMol::AtomName(QString::fromStdString(resinfo->getName()));
                    atomnum = SireMol::AtomNum(resinfo->getSerialNumber());
                    resname = SireMol::ResName(QString::fromStdString(resinfo->getResidueName()));
                    resnum = SireMol::ResNum(resinfo->getResidueNumber());
                    chainname = SireMol::ChainName(QString::fromStdString(resinfo->getChainId()));
                }
                else
                {
                    atomnum = SireMol::AtomNum(n);
                }
            }
            else
            {
                atomnum = SireMol::AtomNum(n);
            }

            // check the molFileAlias property if the atom name hasn't been set
            if (atomname.value().isEmpty())
            {
                std::string alias;

                if (atom->getPropIfPresent<std::string>("molFileAlias", alias))
                {
                    atomname = SireMol::AtomName(QString::fromStdString(alias));
                }
            }

            if (atomname.value().isEmpty())
            {
                atomname = SireMol::AtomName(QString("%1%2").arg(QString::fromStdString(atom->getSymbol())).arg(n));
            }

            // place all atoms into the same cutgroup
            auto a = cg.add(atomnum);

            // now find a residue
            if (resname.value().isEmpty())
            {
                resname = SireMol::ResName("LIG");
            }

            if (resnum.isNull())
            {
                resnum = SireMol::ResNum(1);
            }

            // find this residue - if it doesn't exist, then create it
            SireMol::ResStructureEditor res;

            try
            {
                res = cg.molecule().select(resnum);
            }
            catch (...)
            {
                res = cg.molecule().add(resnum);
                res.rename(resname);
            }

            if (not chainname.value().isEmpty())
            {
                SireMol::ChainStructureEditor chain;

                try
                {
                    chain = cg.molecule().select(chainname);
                }
                catch (...)
                {
                    chain = cg.molecule().add(chainname);
                }

                res.reparent(chain.name());
            }

            a.reparent(res.number());
            a.rename(atomname);

            set_prop(a, "element", SireMol::Element(atom->getAtomicNum()), map);
            set_prop(a, "formal_charge", atom->getFormalCharge() * SireUnits::mod_electron, map);
            set_prop(a, "is_aromatic", qint64(atom->getIsAromatic()), map);
            set_prop(a, "mass", atom->getMass() * SireUnits::g_per_mol, map);
            set_prop(a, "isotope", qint64(atom->getIsotope()), map);
            set_prop(a, "chirality", SireMol::Chirality::fromRDKit(chiral_to_string(atom->getChiralTag())), map);
            set_prop(a, "hybridization", SireMol::Hybridization::fromRDKit(hybridization_to_string(atom->getHybridization())), map);
        }

        // we've built the main structure now - just need to add other
        // properties
        auto molecule = cg.molecule().commit().edit();

        if (map["connectivity"].hasSource())
        {
            auto connectivity = SireMol::Connectivity(molecule.info()).edit();

            const auto bondorder = map["order"];
            const auto stereochemistry = map["stereochemistry"];

            for (const auto &bond : mol->bonds())
            {
                const auto atom0 = atom_to_idx[bond->getBeginAtom()];
                const auto atom1 = atom_to_idx[bond->getEndAtom()];

                connectivity.connect(atom0, atom1);

                const auto b = SireMol::BondID(atom0, atom1);

                if (bondorder.hasSource())
                    connectivity.setProperty(b, bondorder.source(),
                                             SireMol::BondOrder::fromRDKit(bondtype_to_string(bond->getBondType())));

                if (stereochemistry.hasSource())
                    connectivity.setProperty(b, stereochemistry.source(),
                                             SireMol::Stereochemistry::fromRDKit(stereo_to_string(bond->getStereo())));
            }

            molecule.setProperty(map["connectivity"], connectivity.commit());
        }

        if (map["coordinates"].hasSource())
        {
            const int nconformers = mol->getNumConformers();

            if (nconformers > 0)
            {
                const auto &conformer = mol->getConformer(0);

                SireMol::AtomCoords coords(molecule.info());

                for (int i = 0; i < molecule.nAtoms(); ++i)
                {
                    const auto pos = conformer.getAtomPos(i);
                    coords.set(molecule.info().cgAtomIdx(SireMol::AtomIdx(i)),
                               SireMaths::Vector(pos.x, pos.y, pos.z));
                }

                molecule.setProperty(map["coordinates"].source(), coords);
            }
        }

        return molecule.commit();
    }

    SelectorMol rdkit_to_sire(const QList<ROMOL_SPTR> &mols, const PropertyMap &map)
    {
        QList<Molecule> sire_mols;

        for (const auto &mol : mols)
        {
            sire_mols.append(rdkit_to_sire(mol, map));
        }

        return sire_mols;
    }

    QList<ROMOL_SPTR> sire_to_rdkit(const SelectorMol &mols, const PropertyMap &map)
    {
        const int nmols = mols.count();

        if (nmols == 0)
            return QList<ROMOL_SPTR>();

        QVector<ROMOL_SPTR> rdkit_mols(mols.count());

        ROMOL_SPTR *rdkit_mols_data = rdkit_mols.data();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                              {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                rdkit_mols_data[i] = sire_to_rdkit(mols[i], map);
            } });
        }
        else
        {
            for (int i = 0; i < nmols; ++i)
            {
                rdkit_mols_data[i] = sire_to_rdkit(mols[i], map);
            }
        }

        return QList<ROMOL_SPTR>(rdkit_mols.constBegin(), rdkit_mols.constEnd());
    }

    ROMOL_SPTR smarts_to_rdkit(const QString &smarts,
                               const QString &label,
                               const PropertyMap &map)
    {
        RDKit::SmartsParserParams params;
        params.debugParse = 0;
        params.allowCXSMILES = true;
        params.mergeHs = false;
        params.parseName = true;
        params.skipCleanup = false;

        RWMOL_SPTR rdkit_mol;

        if (map.specified("merge_hydrogens"))
        {
            params.mergeHs = map["merge_hydrogens"].value().asABoolean();
        }

        if (map.specified("skip_cleanup"))
        {
            params.skipCleanup = map["skip_cleanup"].value().asABoolean();
        }

        if (map.specified("allow_cx_smiles"))
        {
            params.allowCXSMILES = map["allow_cx_smiles"].value().asABoolean();
        }

        try
        {
            rdkit_mol.reset(RDKit::SmartsToMol(smarts.toStdString(), params));
        }
        catch (...)
        {
        }

        if (rdkit_mol.get() == 0)
        {
            return ROMOL_SPTR();
        }

        rdkit_mol->setProp<std::string>("_Name", label.toStdString());

        return rdkit_mol;
    }

    QList<ROMOL_SPTR> smarts_to_rdkit(const QStringList &smarts,
                                      const QStringList &labels,
                                      const PropertyMap &map)
    {
        if (smarts.count() != labels.count())
            throw SireError::invalid_arg(QObject::tr(
                                             "The number of smarts strings (%1) must match the "
                                             "number of labels (%2)")
                                             .arg(smarts.count())
                                             .arg(labels.count()),
                                         CODELOC);

        const int n = smarts.count();

        if (n == 0)
            return QList<ROMOL_SPTR>();

        QVector<ROMOL_SPTR> ret(n);
        ROMOL_SPTR *ret_data = ret.data();

        if (use_parallel(n, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, n), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    ret_data[i] = smarts_to_rdkit(smarts[i], labels[i], map);
                } });
        }
        else
        {
            for (int i = 0; i < n; ++i)
            {
                ret_data[i] = smarts_to_rdkit(smarts[i], labels[i], map);
            }
        }

        return QList<ROMOL_SPTR>(ret.constBegin(), ret.constEnd());
    }

    ROMOL_SPTR smiles_to_rdkit(const QString &smiles,
                               const QString &label,
                               const PropertyMap &map)
    {
        RDKit::SmilesParserParams params;
        params.debugParse = 0;
        params.sanitize = false;
        params.removeHs = false;
        params.parseName = false;

        RWMOL_SPTR rdkit_mol;

        bool already_sanitized = false;

        if (map.specified("must_sanitize"))
        {
            if (map["must_sanitize"].value().asABoolean())
            {
                // we will force full sanitization
                params.sanitize = true;
                already_sanitized = true;
            }
        }

        try
        {
            rdkit_mol.reset(RDKit::SmilesToMol(smiles.toStdString(), params));
        }
        catch (...)
        {
        }

        if (rdkit_mol.get() == 0)
        {
            return ROMOL_SPTR();
        }

        if (not already_sanitized)
        {
            // Now try to sanitize the molecule - we will do this repeatedly,
            // removing any processes that fail
            unsigned int failed_op;
            unsigned int sanitize_ops = RDKit::MolOps::SANITIZE_ALL;

            bool all_complete = false;
            int ntries = 0;

            while (not all_complete)
            {
                try
                {
                    RDKit::MolOps::sanitizeMol(*rdkit_mol, failed_op, sanitize_ops);
                    all_complete = true;
                }
                catch (...)
                {
                    sanitize_ops |= failed_op;
                }

                ntries += 1;

                if (ntries > 10)
                    // we've tried enough!
                    break;
            }
        }

        rdkit_mol->setProp<std::string>("_Name", label.toStdString());

        bool add_hydrogens = true;
        bool generate_coordinates = true;

        if (map["add_hydrogens"].hasValue())
            add_hydrogens = map["add_hydrogens"].value().asABoolean();

        if (map["generate_coordinates"].hasValue())
            generate_coordinates = map["generate_coordinates"].value().asABoolean();

        if (not add_hydrogens)
            generate_coordinates = false;

        if (add_hydrogens or generate_coordinates)
        {
            try
            {

                RDKit::MolOps::addHs(*rdkit_mol);
            }
            catch (...)
            {
                // could not add hydrogens, so also can't generate coords
                generate_coordinates = false;
            }
        }

        if (generate_coordinates)
        {
            try
            {
                RDKit::DGeomHelpers::EmbedMolecule(*rdkit_mol);
                RDKit::UFF::UFFOptimizeMolecule(*rdkit_mol);
            }
            catch (...)
            {
                // ignore errors
            }
        }

        return rdkit_mol;
    }

    QList<ROMOL_SPTR> smiles_to_rdkit(const QStringList &smiles,
                                      const QStringList &labels,
                                      const PropertyMap &map)
    {
        if (smiles.count() != labels.count())
            throw SireError::invalid_arg(QObject::tr(
                                             "The number of smiles strings (%1) must match the "
                                             "number of labels (%2)")
                                             .arg(smiles.count())
                                             .arg(labels.count()),
                                         CODELOC);

        const int n = smiles.count();

        if (n == 0)
            return QList<ROMOL_SPTR>();

        QVector<ROMOL_SPTR> ret(n);
        ROMOL_SPTR *ret_data = ret.data();

        if (use_parallel(n, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, n), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    ret_data[i] = smiles_to_rdkit(smiles[i], labels[i], map);
                } });
        }
        else
        {
            for (int i = 0; i < n; ++i)
            {
                ret_data[i] = smiles_to_rdkit(smiles[i], labels[i], map);
            }
        }

        return QList<ROMOL_SPTR>(ret.constBegin(), ret.constEnd());
    }

    QString rdkit_to_smarts(const ROMOL_SPTR &mol,
                            const PropertyMap &map)
    {
        if (mol.get() == 0)
            return QString();

        try
        {
            return QString::fromStdString(RDKit::MolToSmarts(*mol));
        }
        catch (...)
        {
            // ignore errors
            return QString();
        }
    }

    QStringList rdkit_to_smarts(const QList<ROMOL_SPTR> &mols,
                                const PropertyMap &map)
    {
        if (mols.isEmpty())
            return QStringList();

        const int nmols = mols.count();

        QVector<QString> smarts(nmols);
        QString *smarts_data = smarts.data();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                              {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                smarts_data[i] = rdkit_to_smarts(mols[i], map);
            } });
        }
        else
        {
            for (int i = 0; i < mols.count(); ++i)
            {
                smarts_data[i] = rdkit_to_smarts(mols[i], map);
            }
        }

        return QStringList(smarts.constBegin(), smarts.constEnd());
    }

    QString rdkit_to_smiles(const ROMOL_SPTR &mol,
                            const PropertyMap &map)
    {
        if (mol.get() == 0)
            return QString();

        try
        {
            return QString::fromStdString(RDKit::MolToSmiles(*mol));
        }
        catch (...)
        {
            // ignore errors
            return QString();
        }
    }

    QStringList rdkit_to_smiles(const QList<ROMOL_SPTR> &mols,
                                const PropertyMap &map)
    {
        if (mols.isEmpty())
            return QStringList();

        const int nmols = mols.count();

        QVector<QString> smiles(nmols);
        QString *smiles_data = smiles.data();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                              {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                smiles_data[i] = rdkit_to_smiles(mols[i], map);
            } });
        }
        else
        {
            for (int i = 0; i < mols.count(); ++i)
            {
                smiles_data[i] = rdkit_to_smiles(mols[i], map);
            }
        }

        return QStringList(smiles.constBegin(), smiles.constEnd());
    }

    ROMOL_SPTR rdkit_remove_hydrogens(const ROMOL_SPTR &mol,
                                      const PropertyMap &map)
    {
        if (mol.get() == 0)
            return mol;

        // really important not to print warnings as these
        // can deadlock jupyter if this is called in parallel!
        RDKit::MolOps::RemoveHsParameters params;
        params.showWarnings = false;

        try
        {
            return ROMOL_SPTR(RDKit::MolOps::removeHs(*mol, params));
        }
        catch (...)
        {
            // ignore errors
            return mol;
        }
    }

    QList<ROMOL_SPTR> rdkit_remove_hydrogens(const QList<ROMOL_SPTR> &mols,
                                             const PropertyMap &map)
    {
        if (mols.isEmpty())
            return QList<ROMOL_SPTR>();

        const int nmols = mols.count();

        QVector<ROMOL_SPTR> ret(nmols);
        ROMOL_SPTR *ret_data = ret.data();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                              {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                ret_data[i] = rdkit_remove_hydrogens(mols[i], map);
            } });
        }
        else
        {
            for (int i = 0; i < mols.count(); ++i)
            {
                ret_data[i] = rdkit_remove_hydrogens(mols[i], map);
            }
        }

        return QList<ROMOL_SPTR>(ret.constBegin(), ret.constEnd());
    }
} // end of namespace SireRDKit

#include "SireMol/select.h"

namespace SireMol
{
    namespace parser
    {

        ////////
        //////// Implementation of the IDSmartsEngine
        ////////

        /** Internal class used to search using smarts strings */
        class IDSmartsEngine : public SelectEngine
        {
        public:
            IDSmartsEngine(bool search_smarts = true);
            SelectEnginePtr createNew(const QList<QVariant> &args) const;

            ~IDSmartsEngine();

            ObjType objectType() const;

        protected:
            SelectResult select(const SelectResult &mols, const PropertyMap &map) const;

            QString smarts;
            bool search_smarts;
        };

        IDSmartsEngine::IDSmartsEngine(bool s) : search_smarts(s)
        {
        }

        IDSmartsEngine::~IDSmartsEngine()
        {
        }

        SelectEngine::ObjType IDSmartsEngine::objectType() const
        {
            return SelectEngine::VIEW;
        }

        SelectEnginePtr IDSmartsEngine::createNew(const QList<QVariant> &args) const
        {
            if (args.isEmpty())
                throw SireError::program_bug(QObject::tr(
                                                 "Weird arguments? %s")
                                                 .arg(Sire::toString(args)),
                                             CODELOC);

            IDSmartsEngine *ptr = new IDSmartsEngine();
            auto p = makePtr(ptr);

            ptr->smarts = args.at(0).toString();
            ptr->search_smarts = this->search_smarts;

            return p;
        }

        SelectResult IDSmartsEngine::select(const SelectResult &mols, const PropertyMap &map) const
        {
            ROMOL_SPTR search_mol;

            if (search_smarts)
            {
                search_mol = SireRDKit::smarts_to_rdkit(smarts, smarts, map);
            }
            else
            {
                search_mol = SireRDKit::smiles_to_rdkit(smarts, smarts, map);
            }

            if (search_mol.get() == 0)
                throw SireError::invalid_key(QObject::tr(
                                                 "Unrecognised smarts string '%1'")
                                                 .arg(smarts),
                                             CODELOC);

            const int match_natoms = search_mol->getNumAtoms();

            const auto water_mask = is_water(mols, map);

            bool water_doesnt_match = false;

            QList<MolViewPtr> ret;

            for (int i = 0; i < mols.count(); ++i)
            {
                const auto mol = mols[i];

                const int natoms = mol.read().nAtoms();

                if (natoms < match_natoms)
                    continue;

                bool this_is_water = water_mask[i];

                if (this_is_water)
                {
                    if (water_doesnt_match)
                        continue;
                }

                auto rdmol = SireRDKit::sire_to_rdkit(mol.read(), map);

                if (rdmol.count() != 1)
                {
                    if (this_is_water)
                        water_doesnt_match = true;

                    continue;
                }

                std::vector<RDKit::MatchVectType> hits;

                if (RDKit::SubstructMatch(*(rdmol.at(0)), *search_mol, hits))
                {
                    if (hits.size() == 0)
                    {
                        if (this_is_water)
                            water_doesnt_match = true;

                        continue;
                    }

                    QList<QList<qint64>> groups;

                    for (const auto &hit : hits)
                    {
                        QList<qint64> selected_atoms;

                        for (const auto &atom : hit)
                        {
                            selected_atoms.append(atom.second);
                        }

                        groups.append(selected_atoms);
                    }

                    ret.append(MolViewPtr(new AtomMatch(mol.read().atoms(), groups)));
                }
                else if (this_is_water)
                {
                    water_doesnt_match = true;
                }
            }

            return SelectResult(ret);
        }

    }
}

namespace SireRDKit
{
    void register_smarts_search()
    {
        SireMol::parser::SelectEngine::registerEngine("smarts",
                                                      new SireMol::parser::IDSmartsEngine(true));

        SireMol::parser::SelectEngine::registerEngine("smiles",
                                                      new SireMol::parser::IDSmartsEngine(false));
    }
}
