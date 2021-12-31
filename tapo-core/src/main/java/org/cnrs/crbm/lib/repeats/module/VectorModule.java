package org.cnrs.crbm.lib.repeats.module;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by pdoviet on 4/15/2015.
 */
public class VectorModule {


    public static void main(String[] args) throws Exception {

        //String pdbCode = "1fdn".toUpperCase();
        //String pdbCode = "2afg".toLowerCase();
//        String pdbCode = "1bhe".toLowerCase();
//        String pdbChain = "A";

        VectorModule vectorModule = new VectorModule();
        System.out.println(vectorModule.extractFeature("1a50", "A"));
//        vectorModule.noTRAnalysis();
//        System.out.println("XXXXXXXXXXXXXXXXXXXXXXXX");
//        vectorModule.repeatsDBAnalysis();

    }


    public void repeatsDBAnalysis() throws Exception {

        StringBuffer buffer = new StringBuffer();

        ProteinCSVReader csvReader = new ProteinCSVReader();
        String out_dir = "C:\\Users\\pdoviet\\Desktop\\";
        String outputfile = "TR.txt";

        List<RowRepeatDB> rows = csvReader.getRepeatDB("data/RDB-dataset.tab");

        for (RowRepeatDB row : rows) {

            if (row.getAnnlevel().equals("Detailed") && !row.getStrclass().equals("II.2")) {
                System.out.println(row.getPdbCode() + "_" + row.getPdbChain() + "\t1");
                //buffer.append(this.extractFeature(row.getPdbCode(), row.getPdbChain()));

            }
        }

        DataIO.writeToFile(buffer.toString(), out_dir + outputfile);
    }

    public void noTRAnalysis() throws Exception {

        List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");
        StringBuffer buffer = new StringBuffer();
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = rows.size();
        for (String row : rows) {


            String pdbCode = row.substring(0, 4);
            String pdbChain = row.substring(5, 6);
            String strClass = row.split("\t")[1];
            if (strClass.equals("1")) {
                //System.out.println(row);
                buffer.append(this.extractFeature(pdbCode, pdbChain));
            }
        }

        String out_dir = "C:\\Users\\pdoviet\\Desktop\\";
        String outputfile = "NOTR.txt";

        DataIO.writeToFile(buffer.toString(), out_dir + outputfile);
    }



    public List<ProVector> extractVectors(Atom[] atoms, String secondaryStr ) throws StructureException {

        List<ProVector> vectors = VectorShape
                .getVectors(atoms, secondaryStr);
        List<ProVector> refined = new ArrayList<ProVector>();
        refined = this.simplifyStructure(vectors);
        for (int i = 0; i < 3; i++) {
            refined = this.simplifyStructure(refined);
        }

        return refined;
    }


    public List<ProVector> extractVectors(String pdbCode, String pdbChain) throws StructureException {


        // String pdbFile = PdbTools.downloadPDB(pdbCode);
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
        // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
        //
        Atom[] atoms = StructureTools.getAtomCAArray(structure
                .getChainByPDB(pdbChain));
        DSSP dssp = new DSSP(pdbCode);
        // make a filter first
        // String secondaryStr = dssp.filterSS(dssp.getSS(atoms));
        String secondaryStr = dssp
                .filterSS(dssp.getCombinedSS(atoms, pdbChain));

        //System.out.println(secondaryStr);
        //VectorModule vectorModule = new VectorModule();


        //Pro2Vector pro2Vector = new Pro2Vector();
        List<ProVector> vectors = VectorShape
                .getVectors(atoms, secondaryStr);

        List<ProVector> refined = new ArrayList<ProVector>();
        refined = this.simplifyStructure(vectors);
        for (int i = 0; i < 3; i++) {
            refined = this.simplifyStructure(refined);
        }


        /**
         * TESTING...........
         */


        /*
        Atom[] ca1 = atoms;
        Atom[] ca2 =  StructureTools.cloneCAArray(ca1);
        CeSymm ceSymm = new CeSymm();
        AFPChain afpChain = ceSymm.pairAlign(ca1, ca2);
        RotationAxis axis = new RotationAxis(afpChain);

        Atom rotationAxis = axis.getRotationAxis();
        Atom rotationPos = axis.getRotationPos();

        // Project each Atom onto the rotation axis to determine limits
        double min, max;
        min = max = Calc.skalarProduct(rotationAxis,atoms[0]);
        for(int i=1;i<atoms.length;i++) {
            double prod = Calc.skalarProduct(rotationAxis,atoms[i]);
            if(prod<min) min = prod;
            if(prod>max) max = prod;
        }
        double uLen = Calc.skalarProduct(rotationAxis,rotationAxis);// Should be 1, but double check
        min/=uLen;
        max/=uLen;

        // Project the origin onto the axis. If the axis is undefined, use the center of mass
        Atom axialPt;
        if(rotationPos == null) {
            Atom center = Calc.centerOfMass(atoms);

            // Project center onto the axis
            Atom centerOnAxis = Calc.scale(rotationAxis, Calc.skalarProduct(center, rotationAxis));

            // Remainder is projection of origin onto axis
            axialPt = Calc.subtract(center, centerOnAxis);

        } else {
            axialPt = rotationPos;
        }

        // Find end points of the rotation axis to display
        Atom axisMin = (Atom) axialPt.clone();
        Calc.scaleAdd(min, rotationAxis, axisMin);
        Atom axisMax = (Atom) axialPt.clone();
        Calc.scaleAdd(max, rotationAxis, axisMax);



        Line line = new Line(new Vector3D(axisMin.getCoords()), new Vector3D((axisMax.getCoords())));

        List<Atom> tmpatomsList = new ArrayList<Atom>();
        for (ProVector v : refined) {
            //if (v.getType().equals("B") || v.getType().equals("H")) {
            tmpatomsList.add(v.getStart());
//            atomsList.add(v.getEnd());
            //}
        }
        List<Atom> atomsList = new ArrayList<Atom>();

        List<Vector3D> point3Ds = new ArrayList<Vector3D>();

        for (int i = 0; i < tmpatomsList.size() - 1; i++) {
            atomsList.add(tmpatomsList.get(i));
            point3Ds.add(new Vector3D(tmpatomsList.get(i).getCoords()));
            Vector3D v1 = new Vector3D(tmpatomsList.get(i).getCoords());
            Vector3D v2 = new Vector3D(tmpatomsList.get(i + 1).getCoords());
            Vector3D v21 = v1.subtract(v2);
            if (v21.getNorm() == 0) {
                i = i + 1;
            }


        }

        for(Vector3D point : point3Ds){
            System.out.print(NumberFormatUtils.format(line.distance(point))+" ");
        }

        */

        return refined;
    }

    public String extractFeature(String pdbCode, String pdbChain) throws StructureException {

        StringBuffer buffer = new StringBuffer();
        int repeatLeng = 0;

        boolean output = false;

        Pro2Vector pro2Vector = new Pro2Vector();
        List<ProVector> refined = this.extractVectors(pdbCode, pdbChain);
        pro2Vector.saveVectorToFile(refined, pdbCode, pdbChain);

        List<Atom> tmpatomsList = new ArrayList<Atom>();
        for (ProVector v : refined) {
            //if (v.getType().equals("B") || v.getType().equals("H")) {
            tmpatomsList.add(v.getStart());
//            atomsList.add(v.getEnd());
            //}
        }
        List<Atom> atomsList = new ArrayList<Atom>();

        List<Vector3D> point3Ds = new ArrayList<Vector3D>();

        for (int i = 0; i < tmpatomsList.size() - 1; i++) {
            atomsList.add(tmpatomsList.get(i));
            point3Ds.add(new Vector3D(tmpatomsList.get(i).getCoords()));
            Vector3D v1 = new Vector3D(tmpatomsList.get(i).getCoords());
            Vector3D v2 = new Vector3D(tmpatomsList.get(i + 1).getCoords());
            Vector3D v21 = v1.subtract(v2);
            if (v21.getNorm() == 0) {
                i = i + 1;
            }

        }


//        Vector3D startPoint = new Vector3D(atomsList.get(0).getCoords());
//        Vector3D endPoint = new Vector3D(atomsList.get(atomsList.size() - 1).getCoords());
//        Line line = new Line(startPoint, endPoint);
//
//
//        double[] his = new double[point3Ds.size()];
//        int index = 0;
//        for (Vector3D point : point3Ds) {
//            his[index] = line.distance(point);
////            System.out.println(line.distance(point));
//            index++;
//        }
//
//
//        SignalProcess signalProcess = new SignalProcess();
//        int minE = 0;
//        int maxE = 0;
//        try {
//            double period = signalProcess.predictPeriod(his);
//            minE = (int) Math.floor(period / 2);
//            maxE = (int) Math.ceil(period / 2);
//            if(maxE==2){
//                minE=2;
//            }
//            if(minE<2)
//                minE =2;
//        } catch (Exception e) {
//            e.printStackTrace();
//        }

        //System.out.println(minE+"XX"+maxE);
        int minE = 2;
        int maxE = 10;
        this.vectorProfile(refined, minE, maxE);


//        Raphael raphael = new Raphael(atomsList,1, 2,0.5);
//        System.out.println(raphael.getRepeatLength()+":"+raphael.getTotalScore());


        //Atom atomSeed = atomsList.get(3);
//        for (int i = 0; i < atomsList.size()-5; i++) {
//
//            Vector3D v1 = new Vector3D(atomsList.get(i).getCoords());
//            Vector3D v2 = new Vector3D(atomsList.get(i + 1).getCoords());
//            Vector3D v3 = new Vector3D(atomsList.get(i + 2).getCoords());
//            Vector3D v4 = new Vector3D(atomsList.get(i + 3).getCoords());
//            Vector3D v5 = new Vector3D(atomsList.get(i + 4).getCoords());
//
//            Vector3D v31 = v1.subtract(v3);
//            Vector3D v35 = v5.subtract(v3);
//
//            double alpha = Dihdral.calcDihedral(v2, v3, v4, v5);
//            double kappa = Math.toDegrees(Vector3D.angle(v31, v35));
//
//            System.out.println(NumberFormatUtils.format(alpha) + "\t" + NumberFormatUtils.format(kappa));
//            //buffer.append(NumberFormatUtils.format(alpha) + "\t" + NumberFormatUtils.format(kappa) + "\n");
//
//            //System.out.print(NumberFormatUtils.format(Calc.getDistance(atomSeed,atomsList.get(i)))+"\t");
//            //System.out.println(v1.getX() + "\t" + v1.getY() + "\t" + v1.getZ());
//
//        }


        //find center point
        double totalX = 0;
        double totalY = 0;
        double totalZ = 0;
        for (int i = 0; i < atomsList.size(); i++) {
            Atom atomi = atomsList.get(i);
            totalX += atomi.getX();
            totalY += atomi.getY();
            totalZ += atomi.getZ();

        }

//        Vector3D centerPoint = new Vector3D(totalX/atomsList.size(),totalY/atomsList.size(), totalZ/atomsList.size());
//
//
//        StringBuffer bf = new StringBuffer();
//        int startIndex = 0;
//        startIndex = atomsList.size() - 1;
//        //startIndex = startIndex/2;
//        Atom atomi = atomsList.get(startIndex);
//        for (int j = 0; j < atomsList.size(); j++) {
//            Atom atomj = atomsList.get(j);
//            Vector3D point = new Vector3D(atomj.getCoords());
//            bf.append((NumberFormatUtils.format(Vector3D.distance(centerPoint,point)) + " "));
//            //bf.append((NumberFormatUtils.format(atomj.getZ()) + " "));
//        }
//
//        System.out.println(bf.toString());




//        String out_dir = "C:\\Users\\pdoviet\\Desktop\\";
//        String outputfile = out_dir + pdbCode.toLowerCase() + pdbChain
//                + "_struct.pdb";
//
//        PdbTools.writePDB(atomsList, outputfile);
        return buffer.toString();

//        List<ProVector> vectorSS = new ArrayList<ProVector>();
//        for (ProVector v : vectors) {
//            if (v.getType().equals("B")) {
//                vectorSS.add(v);
//            }
//
//            vectorsStored.add(v);
//
//        }
//
//
//        for (int i = 0; i < vectorSS.size() - 1; i++) {
//
//            ProVector vi = vectorSS.get(i);
//            ProVector vNext = vectorSS.get(i + 1);
//
//
//            double angle = vectorModule.calsAngle(vi, vNext);
//            double transDist = vectorModule.calsTranslationDis(vi, vNext);
//
//            System.out.println(NumberFormatUtils.format(angle) + "\t" + NumberFormatUtils.format(transDist));
//
//
//        }


//            ProVector v1 = vectors.get(i);
//            ProVector v2 = vectors.get(i + 1);
//            double angle = vectorModule.calsAngle(v1, v2);
//            double dist = vectorModule.calsDistance(v1, v2);
//            System.out.println(NumberFormatUtils.format(angle) + "\t" + NumberFormatUtils.format(dist));


        //List<ProVector> listVectors = vectorModule.getVectors(atoms, secondaryStr);
//
//        ProVector seed = listVectors.get(0);
//        Vector3D startP = new Vector3D(seed.getStart().getCoords());
//        Vector3D endP = new Vector3D(seed.getEnd().getCoords());
//        Vector3D seedVector = endP.subtract(startP);
//
//
//        for (int i = 1; i < atoms.length; i++) {
//            try {
//                Vector3D atom = new Vector3D(atoms[i].getCoords());
//                Vector3D compare = atom.subtract(new Vector3D(atoms[0].getCoords()));
//                double degree = Math.toDegrees(Vector3D.angle(seedVector, compare));
//
//                System.out.print(NumberFormatUtils.format(degree) + ",");
//            } catch (Exception ex) {
//
//            }
//        }


//        for (int i = 0; i < listVectors.size(); i++) {
//
//            ProVector vector = listVectors.get(i);
//            if (vector.getType().equals("H")) {
//                int start = vector.getPosStart();
//                int end = vector.getPosEnd();
//
//                Vector3D startP = new Vector3D(vector.getStart().getCoords());
//                Vector3D endP = new Vector3D(vector.getEnd().getCoords());
//                Vector3D seedVector = endP.subtract(startP);
//
//                //select the middle atom
//                int mid = (end + start) / 2;
//                Atom midAtom = atoms[mid];
//
//                // backward calculation util the before secondary structure or N-terminal
//
//
//                int backwardLimit = start;
//                if (i > 0)
//                    backwardLimit = listVectors.get(i - 1).getPosEnd() + 1;
//                else
//                    backwardLimit = 0;
//
////                for (int j = forwardSearch; j < mid; j++) {
////                    double dist = Calc.getDistance(midAtom, atoms[j]);
////                    System.out.print(NumberFormatUtils.format(dist) + ",");
////                }
//
//
//                // forward calculation util the next secondary structure or C-terminal
//                int forwardLimit = end;
//                if (i < listVectors.size() - 1)
//                    forwardLimit = listVectors.get(i + 1).getPosStart() - 1;
//                else
//                    forwardLimit = atoms.length;
//
//
//                for (int j = backwardLimit; j < forwardLimit * 3; j++) {
//                    //double dist = Calc.getDistance(midAtom, atoms[j]);
//                    //System.out.print(NumberFormatUtils.format(dist) + ",");
//
//                    try {
//                        Vector3D atom = new Vector3D(midAtom.getCoords());
//                        Vector3D compare = new Vector3D(atoms[j].getCoords()).subtract(atom);
//                        double degree = Math.toDegrees(Vector3D.angle(seedVector, compare));
//
//                        System.out.print(NumberFormatUtils.format(degree) + ",");
//                    } catch (Exception ex) {
//
//                    }
//                }
//
//
//                break;
//            }
//
//        }
    }


    private void vectorProfile(List<ProVector> refinedVectors, int minE, int maxE) {

        int MIN_ELEMENTS = minE;
        int MAX_ELEMENTS = maxE;
        Superimposer superimpose = new Superimposer();
        List<ProVector> vectors = ProVector.toSecondaryVector(refinedVectors);
        //List<ProVector> vectors = refinedVectors;
        double maxScore = Double.MIN_VALUE;
        double maxBridgeScore = Double.MIN_VALUE;
        for (int winsize = MIN_ELEMENTS; winsize <= MAX_ELEMENTS; winsize++) {

            //System.out.println();
            List<List<ProVector>> list = new ArrayList<List<ProVector>>();
            for (int i = 0; i < vectors.size() - winsize + 1; i = i + 1) {
                List<ProVector> sub_list = new ArrayList<ProVector>();
                for (int j = 0; j < winsize; j++) {
                    ProVector v = vectors.get(i + j);
                    sub_list.add(v);
                }
                list.add(sub_list);
            }


            for (int j = 0; j < list.size() - winsize; j++) {
                List<Double> scoresList = new ArrayList<Double>();
                List<ProVector> sub_list_ref = list.get(j);
                String pattern_sub_ref = this.getPatterOfVectors(sub_list_ref);
                // right scan
                for (int i = j + winsize; i < list.size() - winsize; i = i + winsize) {
                    List<ProVector> sub_list = list.get(i);
                    String pattern_sub = this.getPatterOfVectors(sub_list);

                    try {
                        if (pattern_sub.equals(pattern_sub_ref)) {
                            double score = superimpose.compareListVector(sub_list_ref,
                                    sub_list);
                            if (score > 0.4)
                                scoresList.add(score);
                            else
                                break;

                        } else {
                            break;
                        }

                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }


                //  left scan
                for (int i = j - winsize; i > 0; i = i - winsize) {
                    List<ProVector> sub_list = list.get(i);
                    String pattern_sub = this.getPatterOfVectors(sub_list);

                    try {
                        if (pattern_sub.equals(pattern_sub_ref)) {
                            double score = superimpose.compareListVector(sub_list_ref,
                                    sub_list);
                            if (score > 0.4)
                                scoresList.add(score);
                            else
                                break;

                        } else {
                            break;
                        }

                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
                double totalScore = 0.0;

                for (double score : scoresList) {
                    totalScore += score;
                }
                if (scoresList.size() > 0) {
                    //totalScore = (double) totalScore / scoresList.size();
                    maxScore = Math.max(totalScore * winsize, maxScore);
                }

            }
            // process score list

            /**
             * find the regions which might be TRs
             */


//            double bridgeScore = 0.0;
//            double maxVScore = 0.0;
//            int countBridge = 0;
//            for (int i = 0; i < scoresList.size(); i++) {
//                if (scoresList.get(i) > 0) {
//                    countBridge++;
//                    maxVScore = Math.max(scoresList.get(i), maxVScore);
//
//
//                    //maxVScore += scoresList.get(i);
//                }
//
//            }
//            bridgeScore = (double) countBridge / scoresList.size();
//            //maxVScore = maxVScore/scoresList.size();
//            maxScore = Math.max((bridgeScore * maxVScore), maxScore);

//            int start = 0;
//            int end = 0;
//            for (int i = 0; i < scoresList.size(); i++) {
//                if (scoresList.get(i) > 0) {
//
//                    start = i;
//                    while (i < scoresList.size()
//                            && scoresList.get(i) > 0) {
//                        i++;
//                    }
//                    i = end = i - 1;
//                    // save
//                    Double totalScore = 0.0;
//                    double max=0.0;
//                    if (end - start > 0) {
//                        List<Double> scores = scoresList.subList(start, end + 1);
//                        for (Double s : scores) {
//                            totalScore += s;
//                            //max = Math.max(s, max);
//                        }
//
//                        //totalScore = max * scores.size();
//                        totalScore= totalScore*winsize;
//                    }
//
//                    maxScore = Math.max(totalScore, maxScore);
//
//                }
//            }


        }

        if (maxScore > Double.MIN_VALUE)
            System.out.println(maxScore);
        else
            System.out.println(0.0);
    }

    private String getPatterOfVectors(List<ProVector> sub_list) {
        String pattern_sub = "";
        for (ProVector v : sub_list) {
            pattern_sub += v.getType();
        }
        return pattern_sub;
    }

    public List<ProVector> simplifyStructure(List<ProVector> vectors) {


        List<ProVector> tmpVectors = new ArrayList<ProVector>();
        //forward
        for (int i = 0; i < vectors.size(); i++) {
            String type = vectors.get(i).getType();
            ProVector vi = vectors.get(i);

            Atom vEnd = (Atom) vi.getEnd().clone();
            //int vPostEnd = vi.getPosStart();
            int vPostEnd = vi.getPosEnd();
            if (type.equals("H") || type.equals("B")) {
                int j = i + 1;
                while (j < vectors.size()) {
                    ProVector vj = vectors.get(j);
                    double angle = this.calsAngle(vi, vj);
                    if (angle < 60) {
                        vEnd = vj.getEnd();
                        //vi.setPosEnd(vj.getPosEnd());
                        vPostEnd = vj.getPosEnd();
                        j++;
                    } else
                        break;
                }

//                vi.setEnd(vEnd);
                vi.getEnd().setCoords(vEnd.getCoords());
                vi.setPosEnd(vPostEnd);
                tmpVectors.add(vi);
                i = j - 1;
            } else
                tmpVectors.add(vi);
        }

        List<ProVector> vectorsStored = new ArrayList<ProVector>();
        //backward
        for (int i = tmpVectors.size() - 1; i >= 0; i--) {
            String type = tmpVectors.get(i).getType();
            ProVector vi = tmpVectors.get(i);
            Atom vStart = (Atom) vi.getStart().clone();
            int vPostStart = vi.getPosStart();
            if (type.equals("H") || type.equals("B")) {
                int j = i - 1;
                while (j > 0) {
                    ProVector vj = tmpVectors.get(j);
                    double angle = this.calsAngle(vi, vj);
                    if (angle < 60) {
                        // fix bugs
//                        vStart = vj.getEnd();
                        vStart = vj.getStart();
                        //vi.setPosEnd(vj.getPosEnd());
//                        vPostStart = vj.getPosEnd();
                        vPostStart = vj.getPosStart();
                        j--;
                    } else
                        break;
                }
//                vi.setStart(vStart);
                vi.getStart().setCoords(vStart.getCoords());
                vi.setPosStart(vPostStart);
                vectorsStored.add(vi);
                i = j + 1;
            } else
                vectorsStored.add(vi);
        }
        Collections.reverse(vectorsStored);

        List<ProVector> connectVectors = new ArrayList<ProVector>();
        // connect vector
        for (int i = 0; i < vectorsStored.size() - 1; i++) {
            connectVectors.add(vectorsStored.get(i));
            ProVector vConnect = new ProVector();
            vConnect.setStart(vectorsStored.get(i).getEnd());
            vConnect.setEnd(vectorsStored.get(i + 1).getStart());
            vConnect.setPosStart(vectorsStored.get(i).getPosEnd());
            vConnect.setPosEnd(vectorsStored.get(i + 1).getPosStart());
            Vector3D v1 = new Vector3D(vConnect.getStart().getCoords());
            Vector3D v2 = new Vector3D(vConnect.getEnd().getCoords());
            if (v1.subtract(v2).getNorm() != 0) {
                vConnect.setType("L");
                connectVectors.add(vConnect);
            }
        }

        if (vectorsStored.size() > 0)
            connectVectors.add(vectorsStored.get(vectorsStored.size() - 1));
        return connectVectors;

    }

    public double calsAngle(ProVector v1, ProVector v2) {

        try {
            Vector3D startP1 = new Vector3D(v1.getStart().getCoords());
            Vector3D endP1 = new Vector3D(v1.getEnd().getCoords());
            Vector3D vector1 = endP1.subtract(startP1);

            Vector3D startP2 = new Vector3D(v2.getStart().getCoords());
            Vector3D endP2 = new Vector3D(v2.getEnd().getCoords());
            Vector3D vector2 = endP2.subtract(startP2);
            return Math.toDegrees(Vector3D.angle(vector1, vector2));
        } catch (Exception ex) {
            return 0;
        }
    }

    public double calsTranslationDis(ProVector v1, ProVector v2) {

        double mid_X1 = (v1.getStart().getX() + v1.getEnd().getX()) / 2;
        double mid_Y1 = (v1.getStart().getY() + v1.getEnd().getY()) / 2;
        double mid_Z1 = (v1.getStart().getZ() + v1.getEnd().getZ()) / 2;

        Vector3D mid1 = new Vector3D(mid_X1, mid_Y1, mid_Z1);


        double mid_X2 = (v2.getStart().getX() + v2.getEnd().getX()) / 2;
        double mid_Y2 = (v2.getStart().getY() + v2.getEnd().getY()) / 2;
        double mid_Z2 = (v2.getStart().getZ() + v2.getEnd().getZ()) / 2;

        Vector3D mid2 = new Vector3D(mid_X2, mid_Y2, mid_Z2);
        return mid2.distance(mid1);


    }

    public double calsDistance(ProVector v1, ProVector v2) {

        Vector3D startP1 = new Vector3D(v1.getStart().getCoords());
        Vector3D endP1 = new Vector3D(v1.getEnd().getCoords());
        Vector3D vector1 = endP1.subtract(startP1);

        Vector3D startP2 = new Vector3D(v2.getStart().getCoords());
        Vector3D endP2 = new Vector3D(v2.getEnd().getCoords());
        Vector3D vector2 = endP2.subtract(startP2);
        return Vector3D.distance(vector1, vector2);
    }

    @Deprecated
    public List<ProVector> getVectors(Atom[] atoms, String secondaryStr) {
        List<ProVector> vectors = new ArrayList<ProVector>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < secondaryStr.length(); i++) {
            char ch = secondaryStr.charAt(i);
            if (ch == 'B') {
                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'B') {
                    i++;
                }
                i = end = i - 1;
                // save
                ProVector v = new ProVector();
                v.setStart(atoms[start]);
                v.setEnd(atoms[end]);
                v.setPosStart(start);
                v.setPosEnd(end);
                v.setType("B");
                int size = end - start + 1;
                vectors.add(v);

            } else if (ch == 'H') {

                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'H') {
                    i++;
                }
                i = end = i - 1;

                ProVector v = new ProVector();
                v.setStart(atoms[start]);
                v.setEnd(atoms[end]);
                v.setPosStart(start);
                v.setPosEnd(end);
                v.setType("H");
                vectors.add(v);
            }


        }

        return vectors;
    }


}
