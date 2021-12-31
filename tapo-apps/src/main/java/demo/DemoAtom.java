package demo;


import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.cnrs.crbm.lib.conf.Dir;

import java.util.List;


/**
 * Created by pdoviet on 3/23/2015.
 */
public class DemoAtom {
    public static void main(String[] args) {

        try {
            PDBFileReader reader = new PDBFileReader();

            // the path to the local PDB installation

            reader.setPath(Dir.PDB_LOCAL);

            // are all files in one directory, or are the files split,
            // as on the PDB ftp servers?
//            reader.setPdbDirectorySplit(true);

            // should a missing PDB id be fetched automatically from the FTP
            // servers?
            reader.setAutoFetch(true);

            reader.setFetchFileEvenIfObsolete(true);
            // should the ATOM and SEQRES residues be aligned when creating
            // the
            // internal data model?
            // reader.setAlignSeqRes(false);

            // should secondary structure get parsed from the file

            FileParsingParameters params = new FileParsingParameters();
            params.setParseSecStruc(true);
            reader.setFileParsingParameters(params);

            Structure struc
                    = reader.getStructureById("3iz8");

            Chain chain = struc.getChainByPDB("A");
            List<Group> groups = chain.getAtomGroups();

            for (Group g : groups) {

                System.out.println(g);
            }


        } catch (Exception ex) {
            // e.printStackTrace();


        }


    }
}
