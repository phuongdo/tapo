package org.cnrs.crbm.trclassification;

import org.apache.commons.io.FileUtils;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by pdoviet on 7/23/2015.
 * Generate PDBs files into groups
 */
public class TRsGenerateFiles {

    public static void main(String[] args) throws IOException, StructureException {
        new TRsGenerateFiles().saveCluster();
    }


    public void saveCluster() throws IOException, StructureException {
        List<String> rows = DataIO.readLines("data/resultLargeScale/trclass-23July.out");
        System.out.println("cleanning output...");
        String dirOut = Dir.CLUSTER_WORKING_DIR + "/output";
        FileUtils.cleanDirectory(new File(dirOut));

        StringBuilder builder = new StringBuilder();
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = rows.size();

        for (String row : rows) {

            bar.update(process, sizeOfProcess);
            process++;
            // do something here
            String[] cols = row.split("\t");
            try {

                String pdb = cols[0];
                String strClass = cols[1];
                double tmScore = Double.parseDouble(cols[3]);
                if (tmScore > 0.5)
                    this.saveFileToOneClass(pdb, strClass);
                else {
                    this.saveFileToOneClass(pdb, "UNK");

                }


            } catch (Exception ex) {

                //ex.printStackTrace();
            }


        }
    }

    public void saveFileToOneClass(String pdb, String strClass) {
        String pdbCode = pdb.substring(0, 4);
        String pdbChain = pdb.substring(5, 6);

        String fileDir = Dir.CLUSTER_WORKING_DIR + "/output/" + strClass;
        File theDir = new File(fileDir);
        if (!theDir.exists()) {
            // make a directory before save file
            theDir.mkdir();
        }

        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        Chain c = null;
        try {
            c = structure.getChainByPDB(pdbChain);
        } catch (StructureException e) {
            // e.printStackTrace();
        }
        Structure newStruct = new StructureImpl();
        newStruct.addChain(c);

        // save file here
        String outputFile = fileDir + "/" + pdbCode + "_" + pdbChain + ".pdb";
        DataIO.writeToFile(newStruct.toPDB(), outputFile);

    }
}
