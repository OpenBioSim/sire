########################################
#
# CMake file for library:SireMol
#
########################################

# Export symbols when building the library
add_definitions( "-DSIREMOL_BUILD" )

# Third Party dependencies of this module

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireMol
set ( SIREMOL_HEADERS
      amberparameters.h
      angleid.h
      atom.h
      atombeads.h
      atomcharges.h
      atomcoords.h
      atomcutting.h
      atomeditor.h
      atomelements.h
      atomenergies.h
      atomforces.h
      atomid.h
      atomidcombos.h
      atomidentifier.h
      atomidx.h
      atomidxmapping.h
      atommapping.h
      atommasses.h
      atommatch.h
      atommatcher.h
      atommatchers.h
      atomname.h
      atomnum.h
      atompolarisabilities.h
      atomproperty.hpp
      atompropertylist.h
      atomradii.h
      atomradicals.h
      atomselection.h
      atomsin.hpp
      atomvelocities.h
      bead.h
      beadeditor.h
      beadid.h
      beadidx.h
      beading.h
      beadnum.h
      beadproperty.hpp
      beads.h
      bondhunter.h
      bondid.h
      bondorder.h
      calculateaxes.hpp
      cgatomidx.h
      cgeditor.h
      cgid.h
      cgidentifier.h
      cgidx.h
      cgname.h
      cgproperty.hpp
      chain.h
      chaineditor.h
      chainid.h
      chainidentifier.h
      chainidx.h
      chainname.h
      chainproperty.hpp
      chainresid.h
      chargeperturbation.h
      chirality.h
      cmapid.h
      connectivity.h
      core.h
      cutgroup.h
      cuttingfunction.h
      dihedralid.h
      editor.hpp
      element-data.h
      element.h
      elementdb.h
      errors.h
      evaluator.h
      evaluatorm.h
      improperid.h
      iswater.h
      geometryperturbation.h
      getghostparam.hpp
      getrmsd.h
      groupatomids.h
      groupgroupids.h
      hybridization.h
      mgid.h
      mgidentifier.h
      mgidsandmaps.h
      mgidx.h
      mgname.h
      mgnum.h
      molatomid.h
      molecule.h
      moleculecutting.h
      moleculedata.h
      moleculegroup.h
      moleculegroups.h
      moleculeinfo.h
      moleculeinfodata.h
      molecules.h
      moleculeview.h
      moleditor.h
      molgroupworkspace.h
      molid.h
      molidentifier.h
      molidx.h
      molinfo.h
      molname.h
      molnum.h
      molresid.h
      molviewproperty.h
      molwithresid.h
      mover.h
      mover.hpp
      mover_metaid.h
      parser.h
      partialmolecule.h
      perturbation.h
      radical.h
      reseditor.h
      resid.h
      residentifier.h
      residue.h
      residuecutting.h
      residx.h
      resin.hpp
      resname.h
      resnum.h
      resproperty.hpp
      segeditor.h
      segid.h
      segidentifier.h
      segidx.h
      segname.h
      segment.h
      segproperty.hpp
      select.h
      selector.hpp
      selectorm.hpp
      selectormol.h
      specifymol.h
      stereochemistry.h
      structureeditor.h
      trajectory.h
      trajectoryaligner.h
      viewsofmol.h
      volumemap.h
      weightfunction.h
      withatoms.h
      within.h
      withres.h
    )

# Define the sources in SireMol
set ( SIREMOL_SOURCES

      register_siremol.cpp

      amberparameters.cpp
      angleid.cpp
      atom.cpp
      atomcoords.cpp
      atomcutting.cpp
      atomeditor.cpp
      atomid.cpp
      atomidentifier.cpp
      atomidxmapping.cpp
      atommapping.cpp
      atommatch.cpp
      atommatcher.cpp
      atommatchers.cpp
      atomproperty.cpp
      atomselection.cpp
      bead.cpp
      beadeditor.cpp
      beadid.cpp
      beading.cpp
      beadproperty.cpp
      beads.cpp
      bondhunter.cpp
      bondid.cpp
      bondorder.cpp
      cgatomidx.cpp
      cgeditor.cpp
      cgid.cpp
      cgidentifier.cpp
      cgproperty.cpp
      chain.cpp
      chaineditor.cpp
      chainid.cpp
      chainidentifier.cpp
      chainproperty.cpp
      chainresid.cpp
      chargeperturbation.cpp
      chirality.cpp
      cmapid.cpp
      connectivity.cpp
      cutgroup.cpp
      cuttingfunction.cpp
      dihedralid.cpp
      elementdb.cpp
      errors.cpp
      evaluator.cpp
      evaluatorm.cpp
      improperid.cpp
      iswater.cpp
      findmcs.cpp
      geometryperturbation.cpp
      getrmsd.cpp
      groupatomids.cpp
      groupgroupids.cpp
      hybridization.cpp
      mgid.cpp
      mgidentifier.cpp
      mgidsandmaps.cpp
      molatomid.cpp
      molgroupworkspace.cpp
      molecule.cpp
      moleculecutting.cpp
      moleculedata.cpp
      moleculegroup.cpp
      moleculegroups.cpp
      moleculeinfo.cpp
      moleculeinfodata.cpp
      molecules.cpp
      moleculeview.cpp
      moleditor.cpp
      molid.cpp
      molidentifier.cpp
      molresid.cpp
      molwithresid.cpp
      molviewproperty.cpp
      mover.cpp
      parser.cpp
      partialmolecule.cpp
      perturbation.cpp
      radical.cpp
      reseditor.cpp
      resid.cpp
      residentifier.cpp
      residue.cpp
      residuecutting.cpp
      resproperty.cpp
      segeditor.cpp
      segid.cpp
      segidentifier.cpp
      segment.cpp
      segproperty.cpp
      select.cpp
      selectormol.cpp
      specifymol.cpp
      stereochemistry.cpp
      structureeditor.cpp
      trajectory.cpp
      trajectoryaligner.cpp
      viewsofmol.cpp
      volumemap.cpp
      weightfunction.cpp
      withatoms.cpp
      within.cpp
      withres.cpp

      ${SIREMOL_HEADERS}
    )

add_library (SireMol ${SIREMOL_SOURCES})

set_target_properties (SireMol PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries, and GSL
target_link_libraries (SireMol
                       SireVol
                       SireCAS
                       SireID
                       SireUnits
                       SireStream
                       )

# installation
install( TARGETS SireMol EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIREMOL_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireMol )
