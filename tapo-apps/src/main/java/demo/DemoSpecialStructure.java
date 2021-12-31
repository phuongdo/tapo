package demo;


import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;

import java.io.IOException;

/**
 * Created by pdoviet on 7/19/2015.
 */
public class DemoSpecialStructure {

    public static void main(String[] args) throws IOException, StructureException {
        // make sure we download chemical component definitions
        // which is required for correctly representing the chromophore
        FileParsingParameters params = new FileParsingParameters();
        params.setLoadChemCompInfo(true);
        params.setAlignSeqRes(true);
        // now register the parameters in the cache
        AtomCache cache = new AtomCache();
        cache.setFileParsingParams(params);
        StructureIO.setAtomCache(cache);


        // request a GFP protein
        Structure s1 = StructureIO.getStructure("1cgd");

        // and print out the internals
        System.out.println(s1.getPDBHeader().toPDB());
        for (Chain c : s1.getChains()) {
            System.out.println("Chain " + c.getChainID() + " details:");
//            System.out.println("Atom ligands: " + c.getAtomLigands());
            System.out.println(c.getSeqResSequence());
        }

        // chromophore is at PDB residue number 66
//        for (Chain c : s1.getChains()) {
//
//            System.out.println("Chain " + c.getChainID() +
//                    " internal " + c.getInternalChainID() +
//                    " ligands " + c.getAtomLigands().size());
//            System.out.println("         10        20        30        40        50        60");
//            System.out.println("1234567890123456789012345678901234567890123456789012345678901234567890");
//            System.out.println(c.getAtomSequence());
//
//            int pos = 0;
//            for (Group g : c.getAtomGroups()) {
//                pos++;
//                System.out.println(pos + " " + g.getResidueNumber() + " " + g.getPDBName() + " " + g.getType() + " " + g.getChemComp().getOne_letter_code() + " " + g.getChemComp().getType());
//            }
//        }
    }
}
