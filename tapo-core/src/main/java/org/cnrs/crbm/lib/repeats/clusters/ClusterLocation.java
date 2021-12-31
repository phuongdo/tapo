package org.cnrs.crbm.lib.repeats.clusters;

import com.apporiented.algorithm.clustering.*;
import org.cnrs.crbm.lib.repeats.RepeatRegion;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 10/10/2014.
 */
public class ClusterLocation {


    public static void main(String[] args) {

        List<Repeat> repeats = new ArrayList<Repeat>();

        Repeat repeat1 = new Repeat();
        repeat1.getRepeats().add(new RepeatContent(307, 374));
        repeat1.getRepeats().add(new RepeatContent(375, 414));
        repeat1.setScore(0.8);
        repeats.add(repeat1);
        Repeat repeat2 = new Repeat();
        repeat2.getRepeats().add(new RepeatContent(307, 374));
        repeat2.getRepeats().add(new RepeatContent(375, 414));
        repeat2.setScore(0.9);
        repeats.add(repeat2);

        System.out.println(repeat1.equals(repeat2));
//        Repeat repeat3 = new Repeat();
//        repeat3.getRepeats().add(new RepeatContent(0, 15));
//        repeat3.setScore(0.2);
//        repeats.add(repeat3);
//        Repeat repeat4 = new Repeat();
//        repeat4.getRepeats().add(new RepeatContent(50, 100));
//        repeat4.setScore(0.82);
//        repeats.add(repeat4);
//        Repeat repeat5 = new Repeat();
//        repeat5.getRepeats().add(new RepeatContent(60, 110));
//        repeat4.setScore(0.81);
//        repeats.add(repeat5);

        List<ClusterRepeat> clusters = new ClusterLocation().cluster(repeats);
        for (ClusterRepeat c : clusters) {
            System.out.println("-----");
            c.sortScoreDESC();
            System.out.println(c);

        }


    }


    private double[][] getDistances(List<Repeat> repeats) {
        int nrows = repeats.size();
        double[][] distances = new double[nrows][nrows];

        for (int i = 0; i < nrows; i++) {
            distances[i][i] = 0.0;
            for (int j = i + 1; j < nrows; j++) {
                Region r1 = new Region(repeats.get(i).getStart(), repeats.get(i).getEnd());
                Region r2 = new Region(repeats.get(j).getStart(), repeats.get(j).getEnd());
                distances[i][j] = 1 - this.overlap(r1, r2);
                distances[j][i] = distances[i][j];
            }
        }
        return distances;

    }

    private static List<String> clusterStrs = new ArrayList<String>();

    public void traverse(Cluster cluster) {

        // Each child of a tree is a root of
        // its subtree.


        if (!cluster.isLeaf()) {

            double distanceC = cluster.getDistance() == null ? 1.0D : cluster
                    .getDistance().doubleValue();
            if (distanceC < 0.5) {
                clusterStrs.add(cluster.toString().replace("Cluster", "").trim());
            } else {

                for (Cluster child : cluster.getChildren()) {
                    double distance = child.getDistance() == null ? 1.0D : child
                            .getDistance().doubleValue();
                    if (distance < 0.5)
                        clusterStrs.add(child.toString().replace("Cluster", "").trim());
                        //System.out.println(child);

                    else
                        traverse(child);


                }


            }
        } else {
            clusterStrs.add(cluster.toString().replace("Cluster", "").trim());
        }
    }


    public List<ClusterRepeat> cluster2times(List<Repeat> repeats) {
        clusterStrs.clear();
        List<ClusterRepeat> clusters = new ArrayList<ClusterRepeat>();
        double[][] distances = this.getDistances(repeats);
        String[] names = new String[repeats.size()];
        for (int i = 0; i < repeats.size(); i++) {
            names[i] = "" + i;
        }
        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
        Cluster cluster = alg.performClustering(distances, names,
                new AverageLinkageStrategy());

        traverse(cluster);

        for (String clusStr : clusterStrs) {
            ClusterRepeat clus = new ClusterRepeat();
            String[] clrs = clusStr.split("&");
            for (String cl_id : clrs)
                clus.getRepeats().add(repeats.get(Integer.parseInt(cl_id)));

            clusters.add(clus);
        }


        return clusters;
    }


    public List<ClusterRepeat> clusterBaseOnPredictedRegion(List<Repeat> repeats, int L) {
        List<ClusterRepeat> clusters = new ArrayList<ClusterRepeat>();
        List<Region> lstPredRegions = ConcensusRegion.findConcensusRegion(repeats, L);

        int label[] = new int[repeats.size()];
        for (Region predRegion : lstPredRegions) {
            ClusterRepeat cluster = new ClusterRepeat();
            cluster.setConcensusRegion(predRegion);
            for (int i = 0; i < repeats.size(); i++) {
                if (label[i] == 0) {
                    Repeat repeat = repeats.get(i);
                    Region repeatRegion = new Region(repeat.getStart(), repeat.getEnd());
                    double cover = this.overlap(repeatRegion, predRegion);
                    // if repeat region cover 50% concensus region
                    if (cover > 0.5) {
                        label[i] = 1;
                        cluster.getRepeats().add(repeat);
                    }
                }
            }
            if (cluster.getRepeats().size() > 0) {
                cluster.assignRScore();
                cluster.sortScoreDESC();
//                cluster.getRepeats().get(0).setCluster("selected");
                clusters.add(cluster);
            }

        }
        return clusters;

    }

    @Deprecated
    public List<ClusterRepeat> cluster(List<Repeat> repeats) {

        List<ClusterRepeat> clusters = new ArrayList<ClusterRepeat>();
        clusters = this.cluster2times(repeats);
        List<Repeat> listRepeats = new ArrayList<Repeat>();
        for (ClusterRepeat c : clusters) {
            // System.out.println("-----");
            //get top
            int minStart = 100;
            int maxEnd = 0;
            // int maxNoTRs = 0;
            double minRU = 1000;
            for (Repeat repeat : c.getRepeats()) {
                int start = repeat.getStart();
                int end = repeat.getEnd();
                minStart = Math.min(start, minStart);
                maxEnd = Math.max(end, maxEnd);
                //maxNoTRs = Math.max(maxNoTRs, repeat.getRepeats().size());
                minRU = Math.min(minRU, repeat.getAvgLength());
            }
            int size = Math.abs(maxEnd - minStart);
            int expectedNoTRs = (int) (size / minRU);

            for (Repeat repeat : c.getRepeats()) {
                int start = repeat.getStart();
                int end = repeat.getEnd();
                repeat.setRankScore(repeat.getScore() * ((double) (end - start + 1) / size) * ((double) repeat.getRepeats().size() / expectedNoTRs));
                //System.out.println(repeat.getRankScore());
            }

            c.sortScoreDESC();

            //get top
            //listRepeats.add(c.getRepeats().get(0));

        }

        //clusters = this.cluster2times(listRepeats);

        return clusters;
    }


    /**
     * Percentage the cover of query with respect to target
     *
     * @param query
     * @param target
     * @return
     */
    public double cover(Region query, Region target) {

        double percent = 0.0;
        int size1 = query.getEnd() - query.getStart() + 1;
        int size2 = target.getEnd() - target.getStart() + 1;
        int size = size2;

        /**
         * RU1 |-------|
         * RU2            |--------|
         */
        if (query.getEnd() < target.getStart())
            return 0.0;


        /**
         * RU1              |-------|
         * RU2  |--------|
         */
        else if (query.getStart() > target.getEnd())
            return 0.0;

        /**
         * RU1  |-------|
         * RU2      |--------|
         */
        else if (target.getStart() <= query.getEnd() && query.getEnd() <= target.getEnd() && query.getStart() < target.getStart()) {
            int match = query.getEnd() - target.getStart() + 1;
            return (double) match / size;
        }
        /**
         * RU1      |-------|
         * RU2  |--------|
         */
        else if (target.getStart() <= query.getStart() && query.getStart() <= target.getEnd() && query.getEnd() > target.getEnd()) {
            int match = target.getEnd() - query.getStart() + 1;
            return (double) match / size;

        }
        /**
         * RU1      |-------|
         * RU2  |-----------------|
         */
        else if (target.getStart() <= query.getStart() && query.getEnd() <= target.getEnd()) {
            return (double) size1 / size;

        }

        /**
         * RU1  |----------------|
         * RU2     |--------|
         */
        else if (query.getStart() <= target.getStart() && target.getEnd() <= query.getEnd()) {

            return 1;

        }

        return 0.0;


    }

    public double overlap(Region repeat1, Region repeat2) {

        // always by smaller
        double percent = 0.0;
        int size1 = repeat1.getEnd() - repeat1.getStart() + 1;
        int size2 = repeat2.getEnd() - repeat2.getStart() + 1;
        int size = Math.min(size1, size2);

        /**
         * RU1 |-------|
         * RU2            |--------|
         */
        if (repeat1.getEnd() < repeat2.getStart())
            return 0.0;


        /**
         * RU1              |-------|
         * RU2  |--------|
         */
        else if (repeat1.getStart() > repeat2.getEnd())
            return 0.0;

        /**
         * RU1  |-------|
         * RU2      |--------|
         */
        else if (repeat2.getStart() <= repeat1.getEnd() && repeat1.getEnd() <= repeat2.getEnd() && repeat1.getStart() < repeat2.getStart()) {
            int match = repeat1.getEnd() - repeat2.getStart() + 1;
            return (double) match / size;
        }
        /**
         * RU1      |-------|
         * RU2  |--------|
         */
        else if (repeat2.getStart() <= repeat1.getStart() && repeat1.getStart() <= repeat2.getEnd() && repeat1.getEnd() > repeat2.getEnd()) {
            int match = repeat2.getEnd() - repeat1.getStart() + 1;
            return (double) match / size;

        }
        /**
         * RU1      |-------|
         * RU2  |-----------------|
         */
        else if (repeat2.getStart() <= repeat1.getStart() && repeat1.getEnd() <= repeat2.getEnd()) {
            return 1;

        }

        /**
         * RU1  |----------------|
         * RU2     |--------|
         */
        else if (repeat1.getStart() <= repeat2.getStart() && repeat2.getEnd() <= repeat1.getEnd()) {

            return 1;

        }

        return 0.0;


    }


}
