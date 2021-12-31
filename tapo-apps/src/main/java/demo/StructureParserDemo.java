package demo;


import org.biojava.nbio.structure.*;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 5/21/2015.
 */
public class StructureParserDemo {

    public static void main(String[] args) throws StructureException {


        String pdbCode = "4dh2";
        String pdbChain = "B";
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        Chain chain = structure.getChainByPDB(pdbChain);
        List<Group> groups = chain.getAtomGroups(GroupType.HETATM);
        String[] atomNames = {"CA"};

        Atom[] atoms = PdbTools.getAtomCAArray(structure.getChainByPDB(pdbChain));


        List<Atom> thisGroupAtoms = new ArrayList<Atom>();
        for (Group g : groups) {


            // a temp container for the atoms of this group

            // flag to check if this group contains all the requested atoms.
            //boolean thisGroupAllAtoms = true;
            for (int i = 0; i < atomNames.length; i++) {
                String atomName = atomNames[i];
                //System.out.println(g);

                Atom a = g.getAtom(atomName);
                thisGroupAtoms.add(a);
            }


        }
        Atom[] calciumsArrs = (Atom[]) thisGroupAtoms.toArray(new Atom[thisGroupAtoms.size()]);

//        for(Atom atom : calciumsArrs){
//            System.out.println(atom);
//            System.out.println(atom.getCoords());
//        }

        for (Atom atom : atoms) {
            System.out.print(NumberFormatUtils.format(Calc.getDistance(atom, calciumsArrs[0])) + " ");
        }


    }
}
