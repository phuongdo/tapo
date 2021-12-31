package demo;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.symmetry.gui.SymmetryDisplay;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.cnrs.crbm.lib.utils.PdbTools;


/**
 * Created by pdoviet on 5/7/2015.
 */
public class DemoCeSymm {


    public static void main(String[] args) throws Exception {
//        System.setProperty(org.slf4j.impl.SimpleLogger.DEFAULT_LOG_LEVEL_KEY, "TRACE");
//
//        final org.slf4j.Logger log = LoggerFactory.getLogger(App.class);
//
//        log.trace("trace");
//        log.debug("debug");
//        log.info("info");
//        log.warn("warning");
//        log.error("error");


//        System.exit(0);
//         new DemoCeSymm().demo("1dvj", "A");
        //new DemoCeSymm().demo("3jut", "A");
        // new DemoCeSymm().demo("1qre", "A");
//        new DemoCeSymm().demo("1qre", "A");
//        new DemoCeSymm().demo("1lxa", "A");
//        new DemoCeSymm().demo("1xhd", "A");
//        new DemoCeSymm().demo("1ap7", "A");
//        new DemoCeSymm().demo("1cwv", "A");

        //new DemoCeSymm().demo("2i13", "A");
        //new DemoCeSymm().demo("2f6e", "A");
//        new DemoCeSymm().demo("4oci", "A");
//        new DemoCeSymm().demo("1tqj", "C");


//        new DemoCeSymm().demo("1l3e", "B");

//        new DemoCeSymm().demo("1ijy", "B");


//        new DemoCeSymm().demo("1lxa", "A");

//        new DemoCeSymm().demo("1u6d", "X");

        new DemoCeSymm().demo("2gsc", "A");


    }

    public void demo(String pdbCode, String pdbChain) throws Exception {


        // String pdbFile = PdbTools.downloadPDB(pdbCode);
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
        // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
        //        System.out.println(structure.toString());


        Atom[] atoms = StructureTools.getRepresentativeAtomArray(structure.getChainByPDB(pdbChain));
        //Atom[] atomsClone = StructureTools.cloneAtomArray(atoms);
        //Initialize the algorithm

//        for (Atom atom : atoms) {
//            System.out.println(atom);
//        }
        CeSymm ceSymm = new CeSymm();

        //Choose some parameters
        CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
        params.setRefineMethod(CESymmParameters.RefineMethod.SINGLE);
//        params.setSymmetryType(CESymmParameters.SymmetryType.CLOSE);
        params.setOptimization(true);
        params.setMultipleAxes(true);


        //Run the symmetry analysis - alignment as an output
        MultipleAlignment symmetry = ceSymm.analyze(atoms, params);


        System.out.println(symmetry.getScore("AvgTM-score"));
        System.out.println(symmetry.size());


        //Test if the alignment returned was refined with
        boolean refined = SymmetryTools.isRefined(symmetry);


        // System.out.println(symmetry.getAtomArrays());

        //  System.out.println(MultipleAlignmentWriter.toAlignedResidues(symmetry));

//        List<List<Integer>> multiples_tmp = MSAWriter.tapoMsaFormat(symmetry);
//        for (List<Integer> lst : multiples_tmp) {
//            System.out.println(lst);
//
//        }


//        List<Atom[]> lstAtoms = symmetry.getAtomArrays();
//        for (Atom[] atomArrs : lstAtoms) {
//
//            for (Atom a : atomArrs) {
//                System.out.print(a.getGroup().getPDBName() + " ");
//            }
//
//            System.out.println();
//
//        }

//Get the axes of symmetry from the aligner
        SymmetryAxes axes = ceSymm.getSymmetryAxes();

//Display the results in jmol with the SymmetryDisplay
        SymmetryDisplay.display(symmetry, axes);

//Show the point group, if any of the internal symmetry
//        QuatSymmetryResults pg = SymmetryTools.getQuaternarySymmetry(symmetry);
//        System.out.println(pg.getSymmetry());


        //
//        Atom[] ca1 = StructureTools.getAtomCAArray(structure
//                .getChainByPDB(pdbChain));
//        Atom[] ca2 = StructureTools.cloneCAArray(ca1);
//
//        String name = pdbCode;
//        CeSymm ceSymm = new CeSymm();
//        AFPChain afpChain = ceSymm.pairAlign(ca1, ca2);
//        afpChain.setName1(name);
//        afpChain.setName2(name);

//        MutilAlign mutilAlign = new MutilAlign();
//        Map<Integer, Integer> pairs = mutilAlign.toAlignedPairs(afpChain, 0, 0);
//        int[] alg1 = new int[pairs.size()];
//        int[] alg2 = new int[pairs.size()];
//        int index = 0;
//        for (Map.Entry<Integer, Integer> entry : pairs.entrySet()) {
//            System.out.println(entry.getKey() + " : "
//                    + entry.getValue());
//            alg1[index] = entry.getKey();
//            alg2[index] = entry.getValue();
//            index++;
//        }
//
//        System.out.println("XX");
//
//        for (int i = 0; i < alg1.length; i=i+20-1) {
//
//            int start = alg1[i];
//            int end = start + 30 - 1;
//            System.out.println(start+":"+end);
//
//        }

        //System.out.println(AfpChainWriter.toDBSearchResult(afpChain));

//        StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain, ca1, ca2);


//
//        int[][][] optAln = afpChain.getOptAln();
//        int[] blockLen = afpChain.getOptLen();
//        for (int block = 0; block < afpChain.getBlockNum(); block++) {
//
//            System.out.println("#######:" + block);
//            for (int i = 0; i < blockLen[block]; i++) {
//                int posA1 = optAln[block][0][i];
//                int posA2 = optAln[block][1][i];
//                System.out.println(posA1 + ":" + posA2);
//
//
//            }
//        }


        //List var12 = DisplayAFP.getPDBresnum(1, afpChain, ca2);

////
//        RotationAxis axis = new RotationAxis(afpChain);
//        jmol.evalString(axis.getJmolScript(ca1));
//        System.out.println(Math.toDegrees(axis.getAngle()));
//        boolean positiveScrew = Math.signum(axis.getRotationAxis().getX()) == Math.signum(axis.getScrewTranslation().getX());
//        System.out.println(positiveScrew);
//
//
//        Atom rotationAxis = axis.getRotationAxis();
//        Atom rotationPos = axis.getRotationPos();
//        Atom[] atoms = ca1;
//
//        // Project each Atom onto the rotation axis to determine limits
//        double min, max;
//        min = max = Calc.skalarProduct(rotationAxis,atoms[0]);
//        for(int i=1;i<atoms.length;i++) {
//            double prod = Calc.skalarProduct(rotationAxis,atoms[i]);
//            if(prod<min) min = prod;
//            if(prod>max) max = prod;
//        }
//        double uLen = Calc.skalarProduct(rotationAxis,rotationAxis);// Should be 1, but double check
//        min/=uLen;
//        max/=uLen;
//
//        // Project the origin onto the axis. If the axis is undefined, use the center of mass
//        Atom axialPt;
//        if(rotationPos == null) {
//            Atom center = Calc.centerOfMass(atoms);
//
//            // Project center onto the axis
//            Atom centerOnAxis = Calc.scale(rotationAxis, Calc.skalarProduct(center, rotationAxis));
//
//            // Remainder is projection of origin onto axis
//            axialPt = Calc.subtract(center, centerOnAxis);
//
//        } else {
//            axialPt = rotationPos;
//        }
//
//        // Find end points of the rotation axis to display
//        Atom axisMin = (Atom) axialPt.clone();
//        Calc.scaleAdd(min, rotationAxis, axisMin);
//        Atom axisMax = (Atom) axialPt.clone();
//        Calc.scaleAdd(max, rotationAxis, axisMax);
//
//
//
//        Line line = new Line(new Vector3D(axisMin.getCoords()), new Vector3D((axisMax.getCoords())));
//        for(Atom atom: ca1){
//            System.out.print(NumberFormatUtils.format(line.distance(new Vector3D(atom.getCoords())))+" ");
//        }
//
//
//
//
//
//        int symmNr = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
//        System.out.println("Symmetry order of: " + symmNr);
//        System.out.println(afpChain.getTMScore());
    }

}
