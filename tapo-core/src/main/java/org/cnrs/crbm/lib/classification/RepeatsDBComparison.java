package org.cnrs.crbm.lib.classification;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.SuperimposerOutput;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 10/20/2015.
 */
public class RepeatsDBComparison {

    List<RowRepeatDB> rowsRef = null;
    List<String> rossList = new ArrayList<String>();

    public static void main(String[] args) {

    }


    public RepeatsDBComparison() {
        try {

            rowsRef = this.getRepeatDB();
            rossList = DataIO.readLines("data/rossmann/rossmannRegions.in");

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Get repeatsDB details annotation.
     *
     * @return
     * @throws Exception
     */


    List<RowRepeatDB> getRepeatDB() throws Exception {
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rows = csvReader.getRepeatDB("data/RDB-dataset.tab");
        List<RowRepeatDB> refs = new ArrayList<RowRepeatDB>();
        for (RowRepeatDB row : rows) {
            if (row.getAnnlevel().equals("Detailed")) {
                refs.add(row);
            }
        }

        return refs;
        // return refs;
    }

    MutilAlign superimpose = new MutilAlign();


    public String classifyProteinBySimilarity(String pdbCode, String pdbChain, String regionAb)
            throws StructureException {

        StringBuffer strBuffer = new StringBuffer();
        int startRegion = Integer.parseInt(regionAb.split("_")[0]);
        int endRegion = Integer.parseInt(regionAb.split("_")[1]);
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        Atom[] atomSet1 = Fragement.getFragementsofAtoms(PdbTools.getAtomCAArray(structure
                .getChainByPDB(pdbChain)), startRegion, endRegion);
        // System.exit(1);
        String assignedBestClass = "UNK";
        double maxTMScore = -1;
        String assignedPdb = "UNK";

        for (RowRepeatDB row : this.rowsRef) {

            try {


//            if(pdb.equals("3m7h_A"))
//                System.out.println();

                String region = row.getRegion();
                String strClass = row.getStrclass();

                int startRegionRef = Integer.parseInt(region.split("-")[0]);
                int endRegionRef = Integer.parseInt(region.split("-")[1]);
                String pdbCodeDB = row.getPdbCode();
                String pdbChainDB = row.getPdbChain();
                Structure structureRef = PdbTools
                        .getStructureFromLocalPdb(pdbCodeDB);

                Atom[] atoms = PdbTools.getAtomCAArray(structureRef
                        .getChainByPDB(pdbChainDB));

                startRegionRef = PdbTools.getPosition(atoms, startRegionRef);
                endRegionRef = PdbTools.getPosition(atoms, endRegionRef);
                Atom[] atomSetRef = Fragement.getFragementsofAtoms(atoms, startRegionRef, endRegionRef);
                // double rmsd = 0.0;
                int len1 = atomSet1.length;
                int len2 = atomSetRef.length;
                // double sim1 = 0.0;
                int avg = (len1 + len2) / 2;
                AFPChain output = superimpose.pairAlignTMScoreDefault(atomSet1, atomSetRef);
                if (output.getTMScore() > maxTMScore) {
                    maxTMScore = output.getTMScore();
                    assignedBestClass = strClass;
                    assignedPdb = row.getPdbCode() + "_" + row.getPdbChain() + "." + region;
                }


            } catch (Exception e) {
                // TODO: handle exception
            }
        }


        // continue to check with rossmann fold

        for (String row : this.rossList) {

            try {
                String[] data = row.split("\t");
                String pdb = data[0];

//            if(pdb.equals("3m7h_A"))
//                System.out.println();
                String strClass = data[2];
                String pdbCodeDB = pdb.substring(0, 4);
                String pdbChainDB = pdb.substring(5, 6);
                Structure structureRef = PdbTools
                        .getStructureFromLocalPdb(pdbCodeDB);

                Atom[] atomSetRef = PdbTools.getAtomCAArray(structureRef
                        .getChainByPDB(pdbChainDB));
                String[] regions = data[1].split(";");
                for (String region : regions) {
                    int startRegionRef = Integer.parseInt(region.split("-")[0]);
                    int endRegionRef = Integer.parseInt(region.split("-")[1]);
                    startRegionRef = PdbTools.getPosition(atomSetRef, startRegionRef);
                    endRegionRef = PdbTools.getPosition(atomSetRef, endRegionRef);
                    Atom[] atomSetRefRegion = Fragement.getFragementsofAtoms(atomSetRef, startRegionRef, endRegionRef);


                    // double rmsd = 0.0;
                    int len1 = atomSet1.length;
                    int len2 = atomSetRef.length;
                    // double sim1 = 0.0;
                    int avg = (len1 + len2) / 2;
                    AFPChain output = superimpose.pairAlignTMScoreDefault(atomSet1, atomSetRefRegion);
                    if (output.getTMScore() > maxTMScore) {
                        maxTMScore = output.getTMScore();
                        assignedBestClass = strClass;
                        //assignedPdb = pdb;
                        assignedPdb = pdbCodeDB + "_" + pdbChainDB + "." + region;
                    }
                }

            } catch (Exception e) {
                // TODO: handle exception
            }
        }
        strBuffer.append(assignedBestClass + "\t" + assignedPdb + "\t" + NumberFormatUtils.format(maxTMScore));
        return strBuffer.toString();
    }

}
