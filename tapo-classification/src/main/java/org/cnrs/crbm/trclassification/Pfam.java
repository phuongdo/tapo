package org.cnrs.crbm.trclassification;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.cath.CathDomain;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.trsfinder.Region;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 10/20/2015.
 */
public class Pfam {


    public static void main(String[] args) {


//         String resNum = "1";
//        Character iCode = 'H';
//        if (iCode == ' ')
//            iCode = null;
//        ResidueNumber residueNumber = new ResidueNumber("A", Integer.valueOf(resNum), iCode);
//        System.out.println(residueNumber.getSeqNum());
//        System.out.println(residueNumber);
//        RepeatFinder repeatFinder = new RepeatFinder("2qrt", "A");
//
//        Atom[] atoms = repeatFinder.getAtoms();
//        for (Atom atom : atoms) {
//
//            System.out.println(atom.getGroup().getResidueNumber().getSeqNum());
//        }

        List<String> rows = DataIO.readLines("input/test.in");
        for (String row : rows) {
            String pdb = row.split("\t")[0];
            String region = row.split("\t")[1];

            String pdbID = pdb.split("\\.")[0];

            String pdbCode = pdbID.substring(0, 4);
            String pdbChain = pdbID.substring(5, 6);

            System.out.println(pdb + "\t" + Pfam.getInstance().assignPfamDomain(pdbCode, pdbChain, region));


        }


    }

    ClusterLocation clusterLocation = new ClusterLocation();

    public String assignPfamDomain(String pdbCode, String pdbChain, String region) {
        // String pfamId = "unk";
        String bestPfamId = "unk";

        double maxOverlap = 0;
        String key = pdbCode + "_" + pdbChain;

        int start;
        int end;
        if (region.startsWith("-")) {
            region = region.substring(1, region.length());
            start = (-1) * Integer.parseInt(region.split("-")[0]);
            end = Integer.parseInt(region.split("-")[1]);

        } else {
            start = Integer.parseInt(region.split("-")[0]);
            end = Integer.parseInt(region.split("-")[1]);
        }


        Region predRegion = new Region(start, end);
        if (pdb_pfam_mapping.containsKey(key)) {
            List<String> rows = pdb_pfam_mapping.get(key);
            for (String row : rows) {
                String[] data = row.split("\t");
                int startRef = Integer.parseInt(data[2].replaceAll("[^0-9]", ""));
                int endRef = Integer.parseInt(data[3].replaceAll("[^0-9]", ""));
                Region refRegion = new Region(startRef, endRef);
//                    trsReferences.getRepeats().add(new RepeatContent(startRef, endRef));
//
                double overlap = clusterLocation.overlap(predRegion, refRegion);
                if (overlap >= 0.2) {

                    if (maxOverlap < overlap) {
                        maxOverlap = overlap;
                        bestPfamId = data[4];
                    }
                    // assignedDomains.add(domain);
                }


            }


        }


        return bestPfamId;
    }

    private static Pfam instance = null;

    Map<String, List<String>> pdb_pfam_mapping = new HashMap<String, List<String>>();


    public Pfam() {

        // load Pfam

        List<String> rows = DataIO.readLines("data/pdb_pfam_mapping.txt");
        for (String row : rows) {
            String[] data = row.split("\t");
            String pdbCode = data[0].toLowerCase();
            String pdbChain = data[1];
            String key = pdbCode + "_" + pdbChain;

            if (pdb_pfam_mapping.containsKey(key)) {
                List<String> lstDomains = pdb_pfam_mapping.get(key);
                lstDomains.add(row);
                pdb_pfam_mapping.put(key, lstDomains);
            } else {
                List<String> lstDomains = new ArrayList<String>();
                lstDomains.add(row);
                pdb_pfam_mapping.put(key, lstDomains);
            }

        }


    }

    public static Pfam getInstance() {
        if (instance == null) {
            instance = new Pfam();
        }
        return instance;
    }

}
