package org.cnrs.crbm.lib.io;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.cnrs.crbm.lib.conf.Dir;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileOutputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.zip.GZIPInputStream;

public class LocalPDBReader {

    static Logger logger = LoggerFactory.getLogger(LocalPDBReader.class);

    public static void main(String[] args) {
        LocalPDBReader loader = new LocalPDBReader();
//        Structure s = loader.getPdb("4hhb");
        Structure s = loader.loadStructureFromFile("input/AF-Q9NQZ3-F1-model_v2.pdb");
//        Structure s = loader.loadStructureFromFile("input/1ib2.pdb");

//        Structure s = loader.loadStructureFromFile("input/rosmann_model.pdb");

        for (Chain c : s.getChains()) {

            System.out.println(">" + c.getChainID());
            String chainID = c.getChainID();
            if (c.getChainID().equals(" "))
                System.out.println("non");
            // print biological sequence
            System.out.println(c.getAtomSequence());
        }

    }

    public Structure loadStructureFromFile(String fileDir) {

        try {
            PDBFileReader reader = new PDBFileReader();

            // the path to the local PDB installation

            reader.setPath(Dir.PDB_LOCAL);

            // are all files in one directory, or are the files split,
            // as on the PDB ftp servers?

//            reader.setPdbDirectorySplit(true);

            reader.setObsoleteBehavior(LocalPDBDirectory.ObsoleteBehavior.FETCH_OBSOLETE);
            // should a missing PDB id be fetched automatically from the FTP
            // servers?
            reader.setAutoFetch(true);
//
//            reader.setFetchFileEvenIfObsolete(true);
            // should the ATOM and SEQRES residues be aligned when creating
            // the
            // internal data model?
            // reader.setAlignSeqRes(false);

            // should secondary structure get parsed from the file

            FileParsingParameters params = new FileParsingParameters();

            params.setParseSecStruc(true);
            //Special Cases When Working with Protein Structures
            // which is required for correctly representing the chromophore
            //A chromophore is the part of a molecule responsible for its color.
            // Some proteins, such as GFP contain a chromopohre that consists of three modified residues.
            // BioJava represents this as a single group in terms of atoms,
            // however as three amino acids when creating the amino acid sequences.
            params.setLoadChemCompInfo(true);
            // should the ATOM and SEQRES residues be aligned when creating the internal data model?
//            params.setAlignSeqRes(true);
            reader.setFileParsingParameters(params);

            Structure structure = reader.getStructure(fileDir);

            return structure;

        } catch (Exception ex) {
            // e.printStackTrace();

            logger.error(ex.getMessage() + " with file name: " + fileDir);

        }

        return null;

    }

    public Structure getPdb(String pdbCode) {

        try {
            PDBFileReader reader = new PDBFileReader(Dir.PDB_LOCAL);

            // the path to the local PDB installation
            System.setProperty("PDB_DIR", Dir.PDB_LOCAL);


            // are all files in one directory, or are the files split,
            // as on the PDB ftp servers?
            //reader.setPdbDirectorySplit(true);


            reader.setObsoleteBehavior(LocalPDBDirectory.ObsoleteBehavior.FETCH_OBSOLETE);
            // should a missing PDB id be fetched automatically from the FTP
            // servers?
//            reader.setAutoFetch(true);
//
//            reader.setFetchFileEvenIfObsolete(true);
            // should the ATOM and SEQRES residues be aligned when creating
            // the
            // internal data model?
            // reader.setAlignSeqRes(false);

            // should secondary structure get parsed from the file

            FileParsingParameters params = new FileParsingParameters();

            params.setParseSecStruc(true);
            //Special Cases When Working with Protein Structures
            // which is required for correctly representing the chromophore
            //A chromophore is the part of a molecule responsible for its color.
            // Some proteins, such as GFP contain a chromopohre that consists of three modified residues.
            // BioJava represents this as a single group in terms of atoms,
            // however as three amino acids when creating the amino acid sequences.
            params.setLoadChemCompInfo(true);
            // should the ATOM and SEQRES residues be aligned when creating the internal data model?
//            params.setAlignSeqRes(true);
            reader.setFileParsingParameters(params);


            Structure structure = reader.getStructureById(pdbCode);

            return structure;

        } catch (Exception ex) {
            // e.printStackTrace();

            logger.error(ex.getMessage() + " with pdb " + pdbCode);

        }

        return null;

    }

    private final static String FTP_PDB = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/obsolete/pdb/";

    public String dowloadPDB(String pdbCode) {
        String fileDir = Dir.PDB_LOCAL + "/" + pdbCode.substring(1, 3) + "/"
                + "pdb" + pdbCode + ".ent.gz";
        // check file exist
        if (!new File(fileDir).exists()) {
            // load
            try {

                String dir = Dir.PDB_LOCAL + "/" + pdbCode.substring(1, 3);
                File theDir = new File(dir);
                if (!theDir.exists()) {
                    // make a directory before save file
                    theDir.mkdir();
                }

                String pdbGZ = FTP_PDB + "/" + pdbCode.substring(1, 3) + "/pdb"
                        + pdbCode + ".ent.gz";
                URL url = new URL(pdbGZ);
                URLConnection con = url.openConnection();
                GZIPInputStream gis = new GZIPInputStream(con.getInputStream());
                FileOutputStream fos = new FileOutputStream(fileDir);
                byte[] buffer = new byte[1024];
                int len;
                while ((len = gis.read(buffer)) != -1) {
                    fos.write(buffer, 0, len);
                }
                // close resources

                fos.close();
                gis.close();

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        return fileDir;

    }
}
