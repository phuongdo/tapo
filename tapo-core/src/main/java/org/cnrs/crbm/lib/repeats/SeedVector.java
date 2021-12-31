package org.cnrs.crbm.lib.repeats;

import java.util.*;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.cnrs.crbm.lib.conf.ThresholdConfig;
import org.cnrs.crbm.lib.math.PeakDetector;
import org.cnrs.crbm.lib.math.Peaks;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.module.VectorWapper;

public class SeedVector {

    public static void vectorsMethod(List<ProVector> vectors) throws Exception {

        int winsize = 2;

        for (int i = 0; i < vectors.size() - winsize; i++) {

            List<ProVector> sub_list = new ArrayList<ProVector>();
            for (int j = 0; j < winsize; j++) {
                ProVector v = vectors.get(i + j);
                sub_list.add(v);
            }
            // find all others has the same distance;
            for (int j = i; j < vectors.size() - winsize; j++) {
                List<ProVector> sub_compare = new ArrayList<ProVector>();
                for (int k = 0; k < winsize; k++) {
                    ProVector v = vectors.get(j + k);
                    sub_compare.add(v);
                }

            }

        }

    }

    public static List<List<Integer>> scanSeed(List<ProVector> vectors)
            throws Exception {

        Superimposer superimposer = new Superimposer();
        List<List<Integer>> seeds = new ArrayList<List<Integer>>();

        int winsize = 2;


        for (int i = 0; i < vectors.size() - winsize; i++) {

            List<Double> signals = new ArrayList<Double>();

            List<ProVector> sub_list = new ArrayList<ProVector>();
            for (int j = 0; j < winsize; j++) {
                ProVector v = vectors.get(i + j);
                sub_list.add(v);
            }

            // find all others has the same distance;
            for (int j = 0; j < vectors.size() - winsize; j++) {
                List<ProVector> sub_compare = new ArrayList<ProVector>();
                for (int k = 0; k < winsize; k++) {
                    ProVector vc = vectors.get(j + k);
                    sub_compare.add(vc);
                }


                String pattern_sub = "";
                for (ProVector v : sub_list) {
                    pattern_sub += v.getType();

                }
                String pattern_com = "";
                for (ProVector v : sub_compare) {
                    pattern_com += v.getType();

                }

                if (pattern_com.equals(pattern_sub))
                    signals.add(superimposer.compareListVector(sub_list,
                            sub_compare));
                else
                    signals.add(0.0);
            }

//
//            for (int j = 0; j < signals.size(); j++) {
//
//                System.out.print(NumberFormatUtils.format(signals.get(j)) + " ");
//            }
//            System.out.println();


            /**
             *convert to position of atom
             */

            List<Integer> seedPos = new ArrayList<Integer>();

//            seedPos.add(i);

            double[] signals_tmp = new double[signals.size()];
            for (int j = 0; j < signals.size(); j++) {
                signals_tmp[j] = signals.get(j);
            }
            LinkedList<Integer> lstPeaks = Peaks.findPeaks(signals_tmp, 2, ThresholdConfig.VECTOR_THRES);
            int[] peaks = new int[lstPeaks.size()];
            for (int j = 0; j < peaks.length; j++) {
                peaks[j] = lstPeaks.get(j);
            }

            //int[] peaks = PeakDetector.getPeaks(signals);
            if (peaks.length > 1) {
                for (int k = 0; k < peaks.length; k++) {
                    //System.out.print(peaks[k] + " ");
                    // from ROC Curve 0.78
                    if (signals.get(peaks[k]) >= ThresholdConfig.VECTOR_THRES) {
                        seedPos.add(vectors.get(peaks[k]).getPosStart());
                    }
                }
            }

            if (seedPos.size() >= 2)
                seeds.add(seedPos);


//            DescriptiveStatistics stats = new DescriptiveStatistics();
//
//            double diff = 0.0;
//            for (int pos = 0; pos < seedPos.size() - 1; pos++) {
//                stats.addValue(seedPos.get(pos + 1) - seedPos.get(pos));
//            }
//            double P_avg = stats.getMean();
//            double P_sd = stats.getStandardDeviation();
//            stats = new DescriptiveStatistics();
//            Set<Integer> setPos = new HashSet<Integer>();
//            for (int pos = 0; pos < seedPos.size() - 1; pos++) {
//                double delta = (seedPos.get(pos + 1) - seedPos.get(pos));
//
//                if (P_avg - P_sd / 2 <= delta
//                        && delta <= P_avg + P_sd / 2) {
//                    stats.addValue(delta);
//                    setPos.add(seedPos.get(pos));
//                    setPos.add(seedPos.get(pos + 1));
//                }
//            }
//
//
//            if (stats.getMean() > 0) {
//
//                for (Integer pos : setPos) {
//                    System.out.print(pos + " ");
//                }
//
//            }
//
//
//            System.out.println("*******************");

        }

        return seeds;

    }

    public static void demoscanSeed(List<ProVector> vectors) {

        // find seed only secondary structure

        // Create and vector x;
        // choose random
        // int random = new Random().nextInt(vectors.size());
        // ProVector vr = vectors.get(0);
        // Vector3D vr1 = new Vector3D(vr.getStart().getCoords());
        // Vector3D vr2 = new Vector3D(vr.getEnd().getCoords());
        //
        // Vector3D v_x = vr2.subtract(vr1);
        Vector3D v_x = new Vector3D(1, 1, 1);
        List<VectorWapper> clusterInput = new ArrayList<VectorWapper>();

        int index = 0;
        for (ProVector v : vectors) {

            if (v.getType().equals("H") || v.getType().equals("B")) {

                // Calculate ange
                Vector3D v1 = new Vector3D(v.getStart().getCoords());
                Vector3D v2 = new Vector3D(v.getEnd().getCoords());
                Vector3D v12 = v2.subtract(v1);
                double angle = Math.toDegrees(Vector3D.angle(v_x, v12));
                clusterInput.add(new VectorWapper(v, angle, index));
            }
            index++;

        }// end for of vectors

        String angleTypesDic = "ABCDEFGHIJKLMN";
        KMeansPlusPlusClusterer<VectorWapper> clusterer = new KMeansPlusPlusClusterer<VectorWapper>(
                4, 10000);
        List<CentroidCluster<VectorWapper>> clusterResults = clusterer
                .cluster(clusterInput);

        // output the clusters
        for (int i = 0; i < clusterResults.size(); i++) {

            List<VectorWapper> parralles_vectors = clusterResults.get(i)
                    .getPoints();
            // for (VectorWapper w : clusterResults.get(i).getPoints()) {
            // if (Math.abs(w.getPoint()[0]
            // - clusterResults.get(i).getCenter().getPoint()[0]) < 5)
            // parralles_vectors.add(w);
            //
            // }

            char angleType = angleTypesDic.charAt(i);

            System.out.println("Cluster " + i + " center:"
                    + clusterResults.get(i).getCenter());
            for (VectorWapper wapper : parralles_vectors) {
                System.out.println(wapper.getVector() + " : "
                        + wapper.getPoint()[0]);
                // if (Math.abs(wapper.getPoint()[0]
                // - clusterResults.get(i).getCenter().getPoint()[0]) < 5)

                vectors.get(wapper.getIndex()).setAngleType(angleType + "");
                // else
                // vectors.get(wapper.getIndex()).setAngleType("x");
            }
            // for (int j = 0; j < parralles_vectors.size() - 1; j++) {
            //
            // ProVector v1 = parralles_vectors.get(j).getVector();
            // ProVector v2 = parralles_vectors.get(j + 1).getVector();
            // Vector3D v11 = new Vector3D(v1.getStart().getCoords());
            // Vector3D v12 = new Vector3D(v1.getEnd().getCoords());
            // Vector3D v21 = new Vector3D(v2.getStart().getCoords());
            // Vector3D v22 = new Vector3D(v2.getEnd().getCoords());
            // Line line1 = new Line(v11, v12);
            // Line line2 = new Line(v21, v22);
            // System.out.println(line2.distance(line1));
            //
            // }

            System.out.println();

        }

        // for (java.util.Map.Entry<Integer, List<ProVector>> a_cluster :
        // clusters
        // .entrySet()) {
        // if (a_cluster.getValue().size() > 2) {
        // System.out.println("****");
        // System.out.println("cluster:" + a_cluster.getKey() + " size:"
        // + a_cluster.getValue().size());
        //
        // List<ProVector> parralles_vectors = a_cluster.getValue();
        //
        // for (int i = 0; i < parralles_vectors.size() - 1; i++) {
        //
        // ProVector v1 = parralles_vectors.get(i);
        // ProVector v2 = parralles_vectors.get(i + 1);
        //
        // Vector3D v11 = new Vector3D(v1.getStart().getCoords());
        // Vector3D v12 = new Vector3D(v1.getEnd().getCoords());
        // Vector3D v21 = new Vector3D(v2.getStart().getCoords());
        // Vector3D v22 = new Vector3D(v2.getEnd().getCoords());
        // Line line1 = new Line(v11, v12);
        // Line line2 = new Line(v21, v22);
        // System.out.println(line1.distance(line2));
        //
        // }
        // }
        //
        // }

    }
}
