package org.cnrs.crbm.lib.dssp;


import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.GzipIO;
import org.cnrs.crbm.lib.io.LocalPDBReader;
import org.cnrs.crbm.lib.utils.ExecUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.io.File;
import java.util.List;

/**
 * This class use for scanning all local install PDB and generating DSSP
 * database
 *
 * @author pdoviet
 */
public class DsspLocalDb {

    LocalPDBReader localPDBReader = new LocalPDBReader();

    public static void main(String[] args) {

        DsspLocalDb dsspLocalDb = new DsspLocalDb();
        dsspLocalDb.generateDB(args[0]);
        // dsspLocalDb.insertDB("1lxa");

    }

    public void generateDB(String fileDir) {
        // read all non-redundancy pdb generated by CDHIT

        //ReadFasta readFasta = new ReadFasta();


        // run new version of dssp binary

        ProgressBar bar = new ProgressBar();

        //System.out.println(Dir.FASTA_LOCAL + "\n");

        try {
            //Set<String> pdbs = readFasta.getAllPdbs(fileDir);
            List<String> pdbs = DataIO.readLines(fileDir);
            System.out.println("Process Starts Now!");
            bar.update(0, pdbs.size());
            int i = 1;
            for (String pdb : pdbs) {
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);
                this.insertDB(pdbCode, pdbChain);
                // update a status bar
                bar.update(i, pdbs.size());
                i++;
            }
            System.out.println("Process Completed!");

        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }


    public void insertDB(String pdbCode, String pdbChain) {

        String fileDir = Dir.DSSP_LOCAL + "/" + pdbCode.substring(1, 3).toLowerCase();
        File theDir = new File(fileDir);
        if (!theDir.exists()) {
            // make a directory before save file
            theDir.mkdir();
        }

        // extract file

        //  String pdbFileGz = Dir.PDB_LOCAL + "/" + pdbCode.substring(1, 3) + "/pdb"
        //           + pdbCode + ".ent.gz";
        String tmpOutFilePDB = Dir.TMP_DIR + "/" + pdbCode + pdbChain + ".pdb";

        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);

        Chain c = null;
        try {
            c = structure.getChainByPDB(pdbChain);
        } catch (StructureException e) {
            // e.printStackTrace();
        }
        Structure newStruct = new StructureImpl();
        newStruct.addChain(c);

        DataIO.writeToFile(newStruct.toPDB(), tmpOutFilePDB);

        // save file here after running DSSP programme
        String outputFile = fileDir + "/" + pdbCode + pdbChain + ".dssp";
        this.calBashshell(tmpOutFilePDB, outputFile);

        // delete tmp file

        new File(tmpOutFilePDB).delete();

    }


    public void generateDSSP(String inputFile, String outputFile) {
        this.calBashshell(inputFile, outputFile);
    }


    public void insertDB(String pdb) {

        String fileDir = Dir.DSSP_LOCAL + "/" + pdb.substring(1, 3).toLowerCase();
        File theDir = new File(fileDir);
        if (!theDir.exists()) {
            // make a directory before save file
            theDir.mkdir();
        }

        // extract file

        String pdbFileGz = Dir.PDB_LOCAL + "/" + pdb.substring(1, 3) + "/pdb"
                + pdb + ".ent.gz";
        String tmpOutFilePDB = Dir.TMP_DIR + "/" + pdb + ".pdb";

        GzipIO.extract(pdbFileGz, tmpOutFilePDB);

        // save file here after running DSSP programme
        String outputFile = fileDir + "/" + pdb + ".dssp";
        this.calBashshell(tmpOutFilePDB, outputFile);

        // delete tmp file

        new File(tmpOutFilePDB).delete();

    }

    private void calBashshell(String inputFile, String outputFile) {

        // new version from 2.0.4.2
        String shellCommand = Dir.DSSP_EXECUTABLE_v2 + " -i " + inputFile
                + " -o " + outputFile;
//        String shellCommand = Dir.DSSP_EXECUTABLE + " " + inputFile + " "
//                + outputFile;
        // System.out.println(shellCommand);
        String strCmdOut = ExecUtils.execShellCmdLinux(shellCommand);

    }

}
