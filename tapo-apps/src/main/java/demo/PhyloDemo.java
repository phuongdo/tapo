package demo;

import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Created by pdoviet on 7/15/2015.
 */
public class PhyloDemo {

    public static void main(final String[] args) {
        // Creating a new rooted tree with two external nodes.
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode d1 = new PhylogenyNode();
        final PhylogenyNode d2 = new PhylogenyNode();
        // Setting of distances.
        d1.setDistanceToParent(1.2);
        d2.setDistanceToParent(4.4);

        // Adding species information.
        final Taxonomy t1 = new Taxonomy();
        t1.setScientificName("TRs");
        d1.getNodeData().addTaxonomy(t1);
        final Taxonomy t2 = new Taxonomy();
        t2.setScientificName("No-TRs");
        d2.getNodeData().addTaxonomy(t2);
        // Adding gene names.
//        final Sequence s1 = new Sequence();
//        s1.setName("Bcl-2");
//        d1.getNodeData().addSequence(s1);
//        final Sequence s2 = new Sequence();
//        s2.setName("Bcl-2");
//        d2.getNodeData().addSequence(s2);
        // Root is a speciation.

        d1.getNodeData().setNodeName("1lxaA");
        d2.getNodeData().setNodeName("1a4yB");
        final Event ev = new Event();
        ev.setSpeciations(1);
        ev.setDuplications(0);
        root.getNodeData().setEvent(ev);
        // Putting the tree together.
        root.addAsChild(d1);
        root.addAsChild(d2);
        phy.setRoot(root);
        phy.setRooted(true);
        // Displaying the newly created tree with Archaeopteryx.
        Archaeopteryx.createApplication(phy);
    }
}
