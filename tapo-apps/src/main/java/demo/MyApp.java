package demo;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.fusesource.jansi.AnsiConsole;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.fusesource.jansi.Ansi.*;
import static org.fusesource.jansi.Ansi.Color.*;

// Import log4j classes.
public class MyApp {

    // Define a static logger variable so that it references the
    // Logger instance named "MyApp".
    static Logger logger = LoggerFactory.getLogger(MyApp.class);

    public static void main(String[] args) {


        int[] a = new int[10];
        for (int i = 0; i < a.length; i++) {
            a[i] = i;
        }


        int winsize = 4;
        for (int i = 4; i < a.length; i = i + winsize) {

            System.out.println(a[i] + ":" + a[i + winsize - 1]);
        }


//        logger.trace("the built-in TRACE level");
//        logger.debug("the built-in DEBUG level");
//        logger.info("the built-in INFO level");
//        logger.warn("the built-in WARN level");
//        logger.error("the built-in ERROR level");
//        System.out.println( ansi().eraseScreen().fg(RED).a("Hello").fg(GREEN).a(" World").reset() );
    }

    public void basicLoad(String pdbId) {
        try {

            PDBFileReader reader = new PDBFileReader();

            // the path to the local PDB installation
            reader.setPath("/tmp");

            // are all files in one directory, or are the files split,
            // as on the PDB ftp servers?
//			reader.setPdbDirectorySplit(true);

            // should a missing PDB id be fetched automatically from the FTP servers?


            // configure the parameters of file parsing

            FileParsingParameters params = new FileParsingParameters();

            // should the ATOM and SEQRES residues be aligned when creating the internal data model?
            params.setAlignSeqRes(true);

            // should secondary structure get parsed from the file
            params.setParseSecStruc(false);

            params.setLoadChemCompInfo(true);

            reader.setFileParsingParameters(params);

            Structure structure = reader.getStructureById(pdbId);

            System.out.println(structure);

            for (Chain c : structure.getChains()) {
                System.out.println("Chain " + c.getChainID() + " details:");
                System.out.println("Atom ligands: " + c.getAtomLigands());
                System.out.println(c.getSeqResSequence());
            }


        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}