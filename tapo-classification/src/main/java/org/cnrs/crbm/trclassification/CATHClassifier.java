package org.cnrs.crbm.trclassification;


import net.sf.javaml.classification.Classifier;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import net.sf.javaml.tools.data.FileHandler;
import net.sf.javaml.tools.weka.WekaClassifier;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathFactory;
import org.biojava.nbio.structure.cath.CathInstallation;
import org.cnrs.crbm.lib.classification.ScopeParser;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import weka.classifiers.bayes.NaiveBayes;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 7/21/2015.
 */
public class CATHClassifier {
    public static void main(String[] args) {

        String pdbID = "4ind";
        String pdbChain = "A";
//
//        CathDatabase cath = CathFactory.getCathDatabase("4.1.0");
//        List<CathDomain> domains = CATHClassifier.getInstance().classifyPdb(pdbID, pdbChain, 0, 0);


        System.out.println(CATHClassifier.getInstance().classifyPdb(pdbID, pdbChain, 0, 0));
//
//        try {
//            // show the structure in 3D
//            for (CathDomain domain : domains) {
//                System.out.println(domain.getClassId());
//            }
//        } catch (Exception e) {
//            // TODO Auto-generated catch block
//            e.printStackTrace();
//        }
    }


    public ProteinCATH classifyPdb(String pdbCode, String pdbChain, int startRegion, int endRegion) {

        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        if (startRegion == endRegion && endRegion == 0) {
            startRegion = 0;
            endRegion = repeatFinder.getAtoms().length - 1;
        }
        Atom[] atoms = repeatFinder.getAtoms();
        //List<CathDomain> assignedDomains = new ArrayList<CathDomain>();
        ProteinCATH proteinCATH = new ProteinCATH();

        CathDomain best = null;
        String strClass = "unk";
        double maxOverlap = 0;
        try {
            List<CathDomain> domains = cath.getDomainsForPdb(pdbCode);
            ClusterLocation clusterLocation = new ClusterLocation();
            Region predRegion = new Region(startRegion, endRegion);


            for (CathDomain domain : domains) {
                for (String rangeCATH : domain.getRanges()) {
                    String pdbChainCATH = rangeCATH.split("_")[0];
                    if (!pdbChainCATH.equals(pdbChain))
                        continue;
                    //String[] rs = rangeCATH.split("_")[1].split("-");
                    Repeat trsFound = new Repeat();
                    trsFound.getRepeats().add(new RepeatContent(startRegion, endRegion));
                    //  Repeat trsReferences = new Repeat();
                    int startRef = 0;
                    int endRef = 0;
                    // convert to reference.
                    String region = rangeCATH.split("_")[1];
                    if (region.startsWith("-")) {
                        region = region.substring(1, region.length());
                        startRef = (-1) * Integer.parseInt(region.split("-")[0]);
                        endRef = Integer.parseInt(region.split("-")[1]);

                    } else {
                        startRef = Integer.parseInt(region.split("-")[0]);
                        endRef = Integer.parseInt(region.split("-")[1]);
                    }

                    startRef = this.getPosition(atoms, startRef);
                    endRef = this.getPosition(atoms, endRef);
                    Region refRegion = new Region(startRef, endRef);
//                    trsReferences.getRepeats().add(new RepeatContent(startRef, endRef));
//
                    double overlap = clusterLocation.overlap(predRegion, refRegion);
                    if (overlap >= 0.2 && pdbChainCATH.equals(pdbChain)) {

                        if (maxOverlap < overlap) {
                            maxOverlap = overlap;
                            best = domain;
                        }
                        // assignedDomains.add(domain);
                    }


                    //if(range.length())
                }


            }


        } catch (Exception ex) {
//            ex.printStackTrace();
        }

        if (maxOverlap >= 0.5 && best != null) {
            proteinCATH.setCathId(best.getCATH());
            proteinCATH.setAssignCathId(best.getClassId() + "");

        } else {

            String ssStr = repeatFinder.getStrSS().substring(startRegion, endRegion + 1);
//                if(pdbCode.equals("1wx0"))
//                    System.out.println();
            String ssPattern = VectorShape.getSSPattern(ssStr);
            double perH = ScopeParser.percentOfSecondaryStructureByType(ssStr, 'H');
            double perB = ScopeParser.percentOfSecondaryStructureByType(ssStr, 'B');
            double perHH = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "HH");
            double perHB = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "HB");
            double perBH = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "BH");
            double perBB = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "BB");
            double[] values = new double[]{perH, perB, perHH, perHB, perBH, perBB};
            Instance instance = new DenseInstance(values);
            strClass = (String) classifier.classify(instance);
            proteinCATH.setAssignCathId(strClass);
        }

        return proteinCATH;

    }


    int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {
            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;
        }
        return 0;
    }

    Classifier classifier = null;
    private CathDatabase cath = null;
    private static CATHClassifier instance = null;

    public CATHClassifier() {

        // load CATHdb classifier

        System.out.println("load CATH classifier");

        cath = new CathInstallation(Dir.CATH_LOCAL);
        cath = CathFactory.getCathDatabase("4.1.0");
        try {
            Dataset data = FileHandler.loadDataset(new File("model/cathv410.train.txt"), 0, "\t");
            NaiveBayes naiveBayes = new NaiveBayes();
            classifier = new WekaClassifier(naiveBayes);//NaiveBayesClassifier(true, true, true);
            classifier.buildClassifier(data);
        } catch (Exception ex) {
            //ex.printStackTrace();

        }
    }

    private Object classify(Instance instance) {
        return classifier.classify(instance);

    }


    public static CATHClassifier getInstance() {
        if (instance == null) {
            instance = new CATHClassifier();
        }
        return instance;
    }

}
