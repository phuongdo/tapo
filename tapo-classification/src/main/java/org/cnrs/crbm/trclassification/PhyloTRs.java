package org.cnrs.crbm.trclassification;

import com.apporiented.algorithm.clustering.*;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

import java.io.FileNotFoundException;
import java.sql.Connection;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 10/5/2015.
 */
public class PhyloTRs {

    private final static String CACHE_MATRIX_FILE = "data/phylo/distMatrix.cache";

    static String TREE_NAME = "tree.ab.in";

    public static void main(String[] args) {
        PhyloTRs phyloTRs = new PhyloTRs("data/phylo/" + TREE_NAME);
        try {


            List<String> pdbIds = phyloTRs.getListOrginalPdbs();
            pdbIds = phyloTRs.preProcessData(pdbIds, "cath");
            pdbIds = phyloTRs.preProcessData(pdbIds, "pfam");
//            for (String pdbId : pdbIds) {
//                System.out.println(pdbId + "\t" + phyloTRs.getPdbTRsFreqs().get(pdbId));
//            }
            phyloTRs.buildPhyloTree(pdbIds);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public PhyloTRs(String treeDir) {
        List<String> rows = DataIO.readLines(treeDir);
        for (String row : rows) {
            String pdbId = row.split("\t")[0];
            // H(omology) level
            //String cathId = row.split("\t")[5];
            String cathId = row.split("\t")[7];
            String pfamId = row.split("\t")[8];
            pdbTRs.put(pdbId, cathId + "|" + pfamId);
            pdbTRsFreqs.put(pdbId, 1);
        }
        tallySeq = DataIO.readLines("data/phylo/tallySeq.verify.txt");


    }


    public List<String> getListOrginalPdbs() {
        List<String> pdbIds = new ArrayList<String>();
        for (Map.Entry<String, String> entry : pdbTRs.entrySet()) {
            pdbIds.add(entry.getKey());
        }
        return pdbIds;
    }

    Map<String, String> pdbTRs = new HashMap<String, String>();
    Map<String, Integer> pdbTRsFreqs = new HashMap<String, Integer>();
    Map<String, String> pdbTRsCATH = new HashMap<String, String>();
    Map<String, Taxonomy> pdbTRsTaxonomy = new HashMap<String, Taxonomy>();
    List<String> tallySeq = new ArrayList<String>();

    public Map<String, Integer> getPdbTRsFreqs() {
        return pdbTRsFreqs;
    }

    public List<String> preProcessData(List<String> pdbIdsInput, String option) {
        System.out.println(option + " pre processing...");
        Map<String, List<String>> pdbGroups = new HashMap<String, List<String>>();
        List<String> pdbIdsUNK = new ArrayList<String>();


        for (String pdbId : pdbIdsInput) {
            String groupKey = pdbTRs.get(pdbId);
            if (option.equals("cath")) {
                groupKey = groupKey.split("\\|")[0];

            } else if (option.equals("pfam")) {
                groupKey = groupKey.split("\\|")[1];
            }

            if (groupKey.equals("unk")) {
                pdbIdsUNK.add(pdbId);

            } else {
                if (pdbGroups.containsKey(groupKey)) {
                    List<String> pdbIds = pdbGroups.get(groupKey);
                    pdbIds.add(pdbId);
                    pdbGroups.put(groupKey, pdbIds);

                } else {
                    List<String> pdbIds = new ArrayList<String>();
                    pdbIds.add(pdbId);
                    pdbGroups.put(groupKey, pdbIds);
                }
            }

        }


//        for (Map.Entry<String, String> entry : pdbTRs.entrySet()) {
//            String groupKey = entry.getValue();
//            if (groupKey.equals("unk")) {
//                pdbIdsUNK.add(entry.getKey());
//
//            } else {
//                if (pdbGroups.containsKey(groupKey)) {
//                    List<String> pdbIds = pdbGroups.get(groupKey);
//                    pdbIds.add(entry.getKey());
//                    pdbGroups.put(groupKey, pdbIds);
//
//                } else {
//                    List<String> pdbIds = new ArrayList<String>();
//                    pdbIds.add(entry.getKey());
//                    pdbGroups.put(groupKey, pdbIds);
//                }
//            }
//
//        }
        List<String> lstPdbIds = new ArrayList<String>();
        for (Map.Entry<String, List<String>> entry : pdbGroups.entrySet()) {
//            System.out.println(entry.getKey());
            List<String> pdbIds = entry.getValue();
            String selectedPdbId = "null";
            int maxLen = -1;
            for (String pdbId : pdbIds) {

                String[] id1 = pdbId.split("\\.");
                int len1 = Integer.parseInt(id1[1].split("_")[1]) - Integer.parseInt(id1[1].split("_")[0]);
                if (len1 > maxLen) {
                    maxLen = len1;
                    selectedPdbId = pdbId;
                }

            }

            pdbTRsFreqs.put(selectedPdbId, pdbIds.size());
            lstPdbIds.add(selectedPdbId);
        }

//        System.out.println("UNK");
        for (String pdbId : pdbIdsUNK) {
//            System.out.println("--- " + pdbId);
            lstPdbIds.add(pdbId);
        }

        return lstPdbIds;

    }

    public void buildPhyloTree(List<String> pdbIds) throws FileNotFoundException {

//        List<String> pdbIds = new ArrayList<String>();
//        for (Map.Entry<String, String> entry : pdbTRs.entrySet()) {
//            pdbIds.add(entry.getKey());
//        }
        System.out.println("starting build a tree...");
        List<String> lstPdbIds = clustering(pdbIds);
        for (int i = 0; i < 50; i++) {
            int before = lstPdbIds.size();
            System.out.println("loop " + (i + 1) + " tree size: " + before);
            lstPdbIds = clustering(lstPdbIds);
            int after = lstPdbIds.size();
            if (before == after) break;
        }


        System.out.println("merging with sequence method from TALY [Richard]...");

        lstPdbIds = this.mergeWithTallyTR(lstPdbIds);
        System.out.println("final tree size: " + lstPdbIds.size());
//        while (lstPdbIds.size() - pdbIds.size() > 0) {
//            pdbIds = lstPdbIds;
//            lstPdbIds = clustering(pdbIds);
//        }
        // save tree
        StringBuffer buffer = new StringBuffer();
        for (String pdbId : lstPdbIds) {
            buffer.append(pdbId + "\n");
        }
        DataIO.writeToFile(buffer.toString(), "data/phylo/" + TREE_NAME + ".out");


        displayTree(lstPdbIds);
    }

    private List<String> mergeWithTallyTR(List<String> lstPdbIds) {

        List<String> filter = new ArrayList<String>();

        for (String pdbId : lstPdbIds) {
            if (tallySeq.contains(pdbId)) {
                filter.add(pdbId);
            }
        }

        return filter;
    }


    private void displayTree(List<String> pdbIds) throws FileNotFoundException {
        // build
        List<String> rows = new ArrayList<String>();
        List<String> rowsName = new ArrayList<String>();
        for (String pdbId : pdbIds) {
            if (pdbTRs.containsKey(pdbId)) {
                rows.add(pdbId);
                String name = this.getPdBSeq(pdbId);
                rowsName.add(name + "|" + pdbTRs.get(pdbId) + "|" + pdbTRsFreqs.get(pdbId));
            }
        }


        String[] names = this.toName(rowsName);
        double[][] distances = this.getDistances(rows);
        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
        Cluster cluster = alg.performClustering(distances, names,
                new AverageLinkageStrategy());

        showPhylogeneticTree(cluster);


    }

    private List<String> clustering(List<String> pdbIds) throws FileNotFoundException {
        // build
        List<String> rows = new ArrayList<String>();
        List<String> rowsName = new ArrayList<String>();
        for (String pdbId : pdbIds) {
            if (pdbTRs.containsKey(pdbId)) {
                rows.add(pdbId);
                //String name = this.getPdBSeq(entry.getKey());
                rowsName.add(pdbId + "|" + pdbTRs.get(pdbId));
            }
        }


        String[] names = this.toName(rowsName);
        double[][] distances = this.getDistances(rows);
        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
        Cluster cluster = alg.performClustering(distances, names,
                new AverageLinkageStrategy());
        List<String> lstPdbIds = new ArrayList<String>();
        traverse(cluster, lstPdbIds);
        return lstPdbIds;

    }

    @Deprecated
    private void hierarchicalClustering() throws FileNotFoundException {

        List<String> rows = new ArrayList<String>();
        List<String> rowsName = new ArrayList<String>();

        for (Map.Entry<String, String> entry : pdbTRs.entrySet()) {
            rows.add(entry.getKey());
            String name = this.getPdBSeq(entry.getKey());
            rowsName.add(name + "|" + entry.getValue());
            pdbTRsCATH.put(name, entry.getValue());


//            if (!pdbTRsTaxonomy.containsKey(entry.getValue())) {
//                final Taxonomy taxonomy = new Taxonomy();
//                taxonomy.setScientificName(entry.getValue());
//                pdbTRsTaxonomy.put(entry.getValue(), taxonomy);
//            }
        }

//        final Taxonomy taxonomyUNK = new Taxonomy();
//        taxonomyUNK.setScientificName("unk");
//        pdbTRsTaxonomy.put("unk", taxonomyUNK);

        //System.exit(1);
        String[] names = this.toName(rowsName);
        double[][] distances = this.getDistances(rows);

//        this.printDistanceMatrix(rows, distances);
        //System.exit(1);
        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
        Cluster cluster = alg.performClustering(distances, names,
                new AverageLinkageStrategy());
//        showPhylogeneticTree(cluster);
//		traverse(cluster);
//        traverse(cluster,l);
//        showPhylogeneticTree(cluster);
    }

    void traverse(Cluster node, List<String> lstPdbIds) {
        if (!node.isLeaf() && node.getChildren().size() == 2 && node.getChildren().get(0).isLeaf() && node.getChildren().get(1).isLeaf()) {
            double distance = node.getDistance() == null ? 100.0D : node
                    .getDistance().doubleValue();
            String[] nodeNames1 = node.getChildren().get(0).getName().split("\\|");
            String[] nodeNames2 = node.getChildren().get(1).getName().split("\\|");
            String[] id1 = nodeNames1[0].split("\\.");
            String[] id2 = nodeNames2[0].split("\\.");
            int len1 = Integer.parseInt(id1[1].split("_")[1]) - Integer.parseInt(id1[1].split("_")[0]);
            int len2 = Integer.parseInt(id2[1].split("_")[1]) - Integer.parseInt(id2[1].split("_")[0]);
//            System.out.println(nodeNames1[0]);
//            System.out.println(nodeNames2[0]);
            if (distance < 0.5 || (!nodeNames1[1].equals("unk") && nodeNames1[1].equals(nodeNames2[1])) || (!nodeNames1[2].equals("unk") && nodeNames1[2].equals(nodeNames2[2]))) {
                if (len1 > len2) {
//                    System.out.println(nodeNames1[0]);
                    String pdbId = nodeNames1[0];
                    lstPdbIds.add(pdbId);
                    pdbTRsFreqs.put(pdbId, pdbTRsFreqs.get(pdbId));
                } else {
//                    System.out.println(nodeNames2[0]);
                    String pdbId = nodeNames2[0];
                    lstPdbIds.add(pdbId);
                    pdbTRsFreqs.put(pdbId, pdbTRsFreqs.get(pdbId));
                }
            } else if ((nodeNames1[1].equals("unk") || nodeNames1[2].equals("unk"))) {
//                System.out.println(nodeNames1[0]);
//                System.out.println(nodeNames2[0]);
                lstPdbIds.add(nodeNames1[0]);
                lstPdbIds.add(nodeNames2[0]);
            } else {
                lstPdbIds.add(nodeNames1[0]);
                lstPdbIds.add(nodeNames2[0]);
            }

        } else if (node.isLeaf()) {
            String[] nodeNames1 = node.getName().split("\\|");
            lstPdbIds.add(nodeNames1[0]);
//            System.out.println(nodeNames1[0]);
        } else {
            // Each child of a tree is a root of
            // its subtree.
            for (Cluster child : node.getChildren()) {
                traverse(child, lstPdbIds);
            }
        }


        // if(child.getDistance())
        // if (child.getDistance() < 30) {

//            if (!child.isLeaf() && child.getChildren().size() == 2 && child.getChildren().get(0).isLeaf() && child.getChildren().get(1).isLeaf()) {
//                double distance = child.getDistance() == null ? 100.0D : child
//                        .getDistance().doubleValue();
//                if (distance < 0.5) {
//                } else {
//
//                }
//
//            } else {
//                traverse(child);
//            }


    }

    private String getPdBSeq(String pdb) {
        String pdbID = pdb.split("\\.")[0];
        String regions = pdb.split("\\.")[1];
        int start = Integer.parseInt(regions.split("_")[0]);
        int end = Integer.parseInt(regions.split("_")[1]);
        String pdbCode = pdbID.substring(0, 4);
        String pdbChain = pdbID.substring(5, 6);
        RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
        start = PdbTools.getResSeq(finder.getAtoms(), start);
        end = PdbTools.getResSeq(finder.getAtoms(), end);
//        if (start < 0) start = 0;
        return pdbID + "." + start + "_" + end;


    }


    private Map<String, Double> getCompareCache() {

        Map<String, Double> map = new HashMap<String, Double>();

        List<String> lines = DataIO.readLines(CACHE_MATRIX_FILE);
        for (String line : lines) {
            String[] rows = line.split("\t");
            map.put(rows[0] + "_" + rows[1], Double.parseDouble(rows[2]));

        }

        return map;

    }

    private double[][] getDistances(List<String> rows)
            throws FileNotFoundException {

        Connection conn = null;
        Statement stmt = null;

        Map<String, Double> cache = this.getCompareCache();
        double[][] distances = new double[rows.size()][rows.size()];
        try {
            for (int i = 0; i < rows.size(); i++) {

                // String pdbCode = rows.get(i).substring(0, 4);
                // String pdbChain = rows.get(i).substring(5, 6);
                for (int j = 0; j < rows.size(); j++) {
                    double score = 1;
                    if (i == j) {
                        score = 0;

                    } else {
                        try {
                            // pdbCode = rows.get(j).substring(0, 4);
                            // pdbChain = rows.get(j).substring(5, 6);
//                            String key1 = rows.get(i).substring(0, 6) + "_" + rows.get(j).substring(0, 6);
//                            String keyInverse = rows.get(j).substring(0, 6) + "_" + rows.get(i).substring(0, 6);
                            String key1 = rows.get(i) + "_" + rows.get(j);
                            String keyInverse = rows.get(j) + "_" + rows.get(i);
                            if (cache.containsKey(key1)) {
                                score = cache.get(key1);
                            } else if (cache.containsKey(keyInverse)) {
                                score = cache.get(keyInverse);
                            }

                        } catch (Exception ex) {
                            // System.out.println(pdbCode);
                            ex.printStackTrace();
                        }

                    }

                    // System.out.println(score);
//					if (score < 0.5)
//						score = 0.0;
//					else
//						score = 1;
                    distances[i][j] = 1 - score;

                }

            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }
        // writer.close();
        return distances;
    }

    private String[] toName(List<String> rows) {

        String[] names = new String[rows.size()];
        for (int i = 0; i < rows.size(); i++) {
            names[i] = rows.get(i).split("\t")[0].trim();
        }
        return names;

    }

    private void showPhylogeneticTree(Cluster cluster) {

        Phylogeny phy = new Phylogeny();
        PhylogenyNode root = new PhylogenyNode();
        buildTree(cluster, root);
        phy.setRoot(root);
        phy.setRooted(true);
        for (final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            //if ( colors.containsKey( n.getName() ) ) {
            // n.getBranchData().setBranchColor( new BranchColor( colors.get( n.getName() ) ) );
            // To make colored subtrees thicker:
            n.getBranchData().setBranchWidth(new BranchWidth(4));
            //}
        }

        //root.getBranchData().setBranchWidth(new BranchWidth(4));


        // Displaying the newly created tree with Archaeopteryx.
        Archaeopteryx.createApplication(phy);


    }

    public void buildTree(Cluster cluster, PhylogenyNode node) {

        for (Cluster child : cluster.getChildren()) {
            PhylogenyNode childNode = new PhylogenyNode();
            //childNode.getBranchData().setBranchWidth(new BranchWidth(4));

            if (!child.isLeaf()) {
                childNode.setName(NumberFormatUtils.format(child.getDistance()));
                node.addAsChild(childNode);
                childNode.setDistanceToParent(child.getDistance());
                buildTree(child, childNode);
            } else {
                childNode.setName(child.getName());
                String pdb = child.getName();//.substring(0, 6);
//                String cathId = "unk";
//                if (pdbTRsCATH.containsKey(pdb)) {
//                    cathId = pdbTRsCATH.get(pdb);
//                }
                // Taxonomy t1 = pdbTRsTaxonomy.get(cathId);
//                childNode.getNodeData().addTaxonomy(pdbTRsTaxonomy.get(cathId));
                node.addAsChild(childNode);
            }

        }

    }
}
