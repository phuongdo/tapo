package org.cnrs.crbm.trclassification;

import net.sf.javaml.classification.Classifier;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import net.sf.javaml.tools.data.FileHandler;
import net.sf.javaml.tools.weka.WekaClassifier;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopInstallation;
import org.cnrs.crbm.lib.classification.ScopeParser;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.trclassification.ProteinSCOP;
import weka.classifiers.bayes.NaiveBayes;

import java.io.File;
import java.util.List;

/**
 * Created by pdoviet on 7/12/2015.
 */
public class ScopeClassifierML {

    private static ScopeClassifierML instance = null;

    Classifier classifier = null;
    // download SCOP if required and load into memory
    ScopInstallation scop = new ScopInstallation(Dir.SCOPE_LOCAL);
    ScopeParser scopeParser = new ScopeParser();

    protected ScopeClassifierML() {
        // Exists only to defeat instantiation.

        try {
            Dataset data = FileHandler.loadDataset(new File("model/scopeDB20.train.txt"), 0, "\t");
            NaiveBayes naiveBayes = new NaiveBayes();
            classifier = new WekaClassifier(naiveBayes);//NaiveBayesClassifier(true, true, true);
            classifier.buildClassifier(data);
        } catch (Exception ex) {

            ex.printStackTrace();

        }
//        double[] values = new double[]{0.761, 0, 0, 0, 0, 0};
//        Instance instance = new DenseInstance(values);
//        Object predictedClassValue = classifier.classify(instance);
//        System.out.println(predictedClassValue);

    }


    public ProteinSCOP classifyPdb(String pdbCode, String pdbChain, int startRegion, int endRegion) {

        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        if (startRegion == endRegion && endRegion == 0) {
            startRegion = 0;
            endRegion = repeatFinder.getAtoms().length - 1;
        }
        ProteinSCOP proteinSCOP = new ProteinSCOP();
        Atom[] atoms = repeatFinder.getAtoms();
        String strClass = "unk";
        try {

            Region predRegion = new Region(startRegion, endRegion);

            List<ScopDomain> domains = scop.getDomainsForPDB(pdbCode);
            ClusterLocation clusterLocation = new ClusterLocation();
            for (ScopDomain domain : domains) {
                for (String range : domain.getRanges()) {
                    String[] rs = range.split(":");
                    if (rs.length > 1 && pdbChain.equals(rs[0])) {
//                        Repeat trsFound = new Repeat();
//                        trsFound.getRepeats().add(new RepeatContent(startRegion, endRegion));
//                        Repeat trsReferences = new Repeat();
                        int startRef = 0;
                        int endRef = 0;
                        // convert to reference.
                        String region = rs[1];
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
//                        trsReferences.getRepeats().add(new RepeatContent(startRef, endRef));

                        Region refRegion = new Region(startRef, endRef);

                        double overlap = clusterLocation.overlap(predRegion, refRegion);
                        if (overlap >= 0.5) {
                            strClass = domain.getClassificationId();
                            break;
                        }

                    } else if (rs.length == 1 && pdbChain.equals(rs[0])) {
                        strClass = domain.getClassificationId();
                        break;


                    }
                    //if(range.length())
                }

//            System.out.println(domain);
//            System.out.println(domain.getRanges());
//            System.out.println(domain.getClassificationId());

            }


            if (strClass.equals("unk")) {
                // continue to classify
                String ssStr = repeatFinder.getStrSS().substring(startRegion, endRegion + 1);

//                if(pdbCode.equals("1wx0"))
//                    System.out.println();
                String ssPattern = VectorShape.getSSPattern(ssStr);
                double perH = scopeParser.percentOfSecondaryStructureByType(ssStr, 'H');
                double perB = scopeParser.percentOfSecondaryStructureByType(ssStr, 'B');
                double perHH = scopeParser.percentageArrangemanceOfStructure(ssPattern, "HH");
                double perHB = scopeParser.percentageArrangemanceOfStructure(ssPattern, "HB");
                double perBH = scopeParser.percentageArrangemanceOfStructure(ssPattern, "BH");
                double perBB = scopeParser.percentageArrangemanceOfStructure(ssPattern, "BB");
                double[] values = new double[]{perH, perB, perHH, perHB, perBH, perBB};
                Instance instance = new DenseInstance(values);
//                StringBuffer builder = new StringBuffer();
//                builder.append(pdbCode + ",");
//                builder.append(NumberFormatUtils.format(perH) + ",");
//                builder.append(NumberFormatUtils.format(perB) + ",");
//                builder.append(NumberFormatUtils.format(perHH) + ",");
//                builder.append(NumberFormatUtils.format(perHB) + ",");
//                builder.append(NumberFormatUtils.format(perBH) + ",");
//                builder.append(NumberFormatUtils.format(perBB));
//                System.out.println(builder.toString());
                strClass = (String) classifier.classify(instance);

                proteinSCOP.setAssignScopId(strClass);
            } else {
                proteinSCOP.setScopId(strClass);
                proteinSCOP.setAssignScopId(strClass.split("\\.")[0]);

            }
        } catch (Exception ex) {
            // do something here :)
        }

        return proteinSCOP;
    }


    int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {
            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;
        }
        return 0;
    }

    private Object classify(Instance instance) {
        return classifier.classify(instance);

    }

    public static ScopeClassifierML getInstance() {
        if (instance == null) {
            instance = new ScopeClassifierML();
        }
        return instance;
    }
}
