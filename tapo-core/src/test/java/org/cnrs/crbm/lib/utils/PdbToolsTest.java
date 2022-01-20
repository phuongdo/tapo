package org.cnrs.crbm.lib.utils;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class PdbToolsTest {

    @Test
    public void testLoadStructureFromFile() throws Exception{
        Structure s = PdbTools.getStructureFromFile("/Users/phuongdv/Developer/CNRS/tapo/input/AF-Q9NQZ3-F1-model_v2.pdb");
//        Structure s = PdbTools.getStructureFromFile("/Users/phuongdv/Developer/CNRS/tapo/input/1ib2.pdb");

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

}