package org.cnrs.crbm.trclassification;

import org.cnrs.crbm.lib.io.DataIO;
import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by pdoviet on 7/25/2015.
 */
public class PhyloCATH {

    public static void main(String[] args) throws Exception {
        new PhyloCATH().buildTree();
    }

    public void buildTree() throws Exception {
        // read the file

        String fileDir = "data/trclassification/UnkProteinTRs";
        List<String> lines = DataIO.readLines(fileDir);


        // build a roots

        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();


        // highest classes are length?
        List<PhylogenyNode> lenghNodes = new ArrayList<PhylogenyNode>();


        //
        List<PhylogenyNode> cNodes = new ArrayList<PhylogenyNode>();
        List<PhylogenyNode> aNodes = new ArrayList<PhylogenyNode>();
        List<PhylogenyNode> tNodes = new ArrayList<PhylogenyNode>();
        List<PhylogenyNode> hNodes = new ArrayList<PhylogenyNode>();

        Set<String> setLengthLevels = new HashSet<String>();

        Set<String> setCLevels = new HashSet<String>();
        Set<String> setALevels = new HashSet<String>();
        Set<String> setTLevels = new HashSet<String>();
        Set<String> setHLevels = new HashSet<String>();
        for (String line : lines) {
            String[] cols = line.split("\t");
            String pdb = cols[0];
            String _C = cols[1];
            String _A = cols[2];
            String _T = cols[3];
            String _H = cols[4];
            Double avgL = Double.parseDouble(cols[5]);

            String _L = "1";
            if (avgL <= 10) {
                _L = "10L";
            } else if (avgL <= 20) {
                _L = "20L";
            } else if (avgL <= 30) {
                _L = "30L";
            } else if (avgL <= 40) {
                _L = "40L";
            } else if (avgL <= 50) {
                _L = "50L";
            } else if (avgL <= 60) {
                _L = "60L";
            } else if (avgL <= 70) {
                _L = "70L";
            } else if (avgL <= 80) {
                _L = "80L";
            } else if (avgL <= 90) {
                _L = "90L";
            } else if (avgL > 90) {
                _L = "X90L";
            }


            if (!setLengthLevels.contains(_L)) {
                setLengthLevels.add(_L);
                PhylogenyNode node = new PhylogenyNode();
                node.setName(_L);
                lenghNodes.add(node);
            }

            if (!setCLevels.contains(_L + "." + _C)) {
                setCLevels.add(_L + "." + _C);
                PhylogenyNode node = new PhylogenyNode();
                node.setName(_L + "." + _C);
                cNodes.add(node);
            }
            if (!setALevels.contains(_L + "." + _C + "." + _A)) {
                setALevels.add(_L + "." + _C + "." + _A);
                PhylogenyNode node = new PhylogenyNode();
                node.setName(_L + "." + _C + "." + _A);
                aNodes.add(node);
            }
            if (!setTLevels.contains(_L + "." + _C + "." + _A + "." + _T)) {
                setTLevels.add(_L + "." + _C + "." + _A + "." + _T);
                PhylogenyNode node = new PhylogenyNode();
                node.setName(_L + "." + _C + "." + _A + "." + _T);
                tNodes.add(node);
            }
            if (!setHLevels.contains(_L + "." + _C + "." + _A + "." + _T + "." + _H)) {
                setHLevels.add(_L + "." + _C + "." + _A + "." + _T + "." + _H);
                PhylogenyNode node = new PhylogenyNode();
                node.setName(_L + "." + _C + "." + _A + "." + _T + "." + _H);
                hNodes.add(node);
            }


        }

        // Step 1. Add leaves

        // join to H level

        for (String line : lines) {
            String[] cols = line.split("\t");
            String pdb = cols[0];
            String _C = cols[1];
            String _A = cols[2];
            String _T = cols[3];
            String _H = cols[4];

            Double avgL = Double.parseDouble(cols[5]);

            String _L = "1";
            if (avgL <= 10) {
                _L = "10L";
            } else if (avgL <= 20) {
                _L = "20L";
            } else if (avgL <= 30) {
                _L = "30L";
            } else if (avgL <= 40) {
                _L = "40L";
            } else if (avgL <= 50) {
                _L = "50L";
            } else if (avgL <= 60) {
                _L = "60L";
            } else if (avgL <= 70) {
                _L = "70L";
            } else if (avgL <= 80) {
                _L = "80L";
            } else if (avgL <= 90) {
                _L = "90L";
            } else if (avgL > 90) {
                _L = "X90L";
            }

            for (PhylogenyNode node : hNodes) {

                if (node.getName().equals(_L + "." + _C + "." + _A + "." + _T + "." + _H)) {
                    PhylogenyNode leaf = new PhylogenyNode();
                    leaf.getNodeData().setNodeName(pdb);
                    node.addAsChild(leaf);
                }

            }
        }
        // join to T level
        for (PhylogenyNode tnode : tNodes) {
            for (PhylogenyNode hnode : hNodes) {
                String[] cols = hnode.getName().split("\\.");
                String _L = cols[0];
                String _C = cols[1];
                String _A = cols[2];
                String _T = cols[3];
                if (tnode.getName().equals(_L + "." + _C + "." + _A + "." + _T)) {
                    tnode.addAsChild(hnode);
                }
            }

        }

        // join to A level
        for (PhylogenyNode anode : aNodes) {
            for (PhylogenyNode tnode : tNodes) {
                String[] cols = tnode.getName().split("\\.");
                String _L = cols[0];
                String _C = cols[1];
                String _A = cols[2];

                if (anode.getName().equals(_L + "." + _C + "." + _A)) {
                    anode.addAsChild(tnode);
                }
            }

        }

        // join to C level
        for (PhylogenyNode cnode : cNodes) {
            for (PhylogenyNode anode : aNodes) {
                String[] cols = anode.getName().split("\\.");
                String _L = cols[0];
                String _C = cols[1];
                if (cnode.getName().equals(_L + "." + _C)) {
                    cnode.addAsChild(anode);
                }
            }

        }

        // join to L level
        for (PhylogenyNode lnode : lenghNodes) {
            for (PhylogenyNode cnode : cNodes) {
                String[] cols = cnode.getName().split("\\.");
                String _L = cols[0];
                if (lnode.getName().equals(_L)) {
                    lnode.addAsChild(cnode);
                }
            }
            root.addAsChild(lnode);
        }

        phy.setRoot(root);
        phy.setRooted(true);
        Archaeopteryx.createApplication(phy);
    }


    private void addLeafToTree(PhylogenyNode leaf, PhylogenyNode root) {
        // level H
        for (PhylogenyNode child : root.getDescendants()) {
            child.setName("");
        }

    }
}
