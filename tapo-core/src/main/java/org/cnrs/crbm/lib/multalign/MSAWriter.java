package org.cnrs.crbm.lib.multalign;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by pdoviet on 9/7/2015.
 */
public class MSAWriter {
    public MSAWriter() {
    }


    public static List<List<Integer>> tapoMsaFormat(MultipleAlignment multAln) {

        List<List<Integer>> msa = new ArrayList<List<Integer>>();
        for (int str = 0; str < multAln.size(); ++str) {
            msa.add(new ArrayList<Integer>());
        }
//        int nRow = multAln.size();
//        int nCol = multAln.getCoreLength();
//        int[][] multi = new int[nRow][nCol];

        int i$;


//        for (i$ = 0; i$ < multAln.size(); ++i$) {
//            residueGroup.append("#Num" + (i$ + 1) + "\tChain" + (i$ + 1) + "\tAA" + (i$ + 1) + "\t");
//        }
//
//        residueGroup.append("\n");
        Iterator var8 = multAln.getBlocks().iterator();

        while (var8.hasNext()) {
            Block b = (Block) var8.next();
            for (int res = 0; res < b.length(); ++res) {
                for (int str = 0; str < multAln.size(); ++str) {
//                    int start = Integer.parseInt(units[str].split("-")[0]);

                    Integer residue = (Integer) ((List) b.getAlignRes().get(str)).get(res);

                    int pos = -1;
                    if (residue == null) {
                        pos = -1;
                    } else {
                        //Atom atom = ((Atom[]) multAln.getAtomArrays().get(str))[residue.intValue()];
                        pos = residue.intValue();
                    }

                    msa.get(str).add( pos);
                }

            }
        }

        return msa;
    }

}
