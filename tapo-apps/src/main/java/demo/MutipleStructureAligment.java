package demo;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcMain;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcParameters;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.cnrs.crbm.lib.multalign.MSAWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by pdoviet on 9/3/2015.
 */
public class MutipleStructureAligment {
    public static void main(String[] args) throws StructureException, IOException {

        //Specify the structures to pairAlign: some ASP-proteinases
        List<String> names = Arrays.asList("3app", "4ape", "5pep");

//Load the CA atoms of the structures
        AtomCache cache = new AtomCache();
        List<Atom[]> atomArrays = new ArrayList<Atom[]>();
        for (String name : names) {
            atomArrays.add(cache.getAtoms(name));
        }

        //Here the multiple structural alignment algorithm comes in place to generate the alignment object
        MultipleMcMain algorithm = new MultipleMcMain(new CeMain());
        MultipleMcParameters params = (MultipleMcParameters) algorithm.getParameters();
        params.setMinBlockLen(15);
        params.setMinAlignedStructures(10);

        MultipleAlignment result = algorithm.align(atomArrays);
        result.getEnsemble().setStructureNames(names);

        //Information about the alignment
        result.getEnsemble().setAlgorithmName(algorithm.getAlgorithmName());
        result.getEnsemble().setVersion(algorithm.getVersion());
//        System.out.println(MultipleAlignmentWriter.toAlignedResidues(result.getEnsemble()));
        //Output the sequence alignment + transformations
//        System.out.println(MSAWriter.toFatCat(result));
//        System.out.println(MSAWriter.tapoMsaFormat(result));
//        System.out.println(MultipleAlignmentWriter.toTransformMatrices(result));
//        System.out.println(MultipleAlignmentWriter.toXML(result.getEnsemble()));
//        MultipleAlignmentDisplay.display(result);

    }
}
