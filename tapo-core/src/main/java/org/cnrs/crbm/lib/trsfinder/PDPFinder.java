package org.cnrs.crbm.lib.trsfinder;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.domain.pdp.CutDomain;
import org.biojava.nbio.structure.domain.pdp.Domain;
import org.biojava.nbio.structure.domain.pdp.Segment;
import org.cnrs.crbm.lib.domain.LocalProteinDomainParser;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.SuperimposerOutput;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

public class PDPFinder extends Finder {

    // Atom[] atoms = null;

    static Logger logger = LoggerFactory.getLogger(PDPFinder.class);
    public PDPFinder(Features features) {
        super(features);
        this.name = "RMSD";
        // atoms = features.getAtoms();
    }


    static Superimposer superimposer = new Superimposer();
    double maxTM = 0.0;
    @Override
    public void findRepeat(Features features) {
        try {

            CutDomain.verbose = false;
            List<Domain> domains = LocalProteinDomainParser
                    .suggestDomains(atoms);
            String ss = features.getStrSS();
            domains = this.removeOverlapDomain(domains);
            //
            // System.out.println("RESULTS: =====");

//            for (Domain domain : domains) {
//                List<Segment> segments = domain.getSegments();
//                System.out.println(segments.get(0).getFrom()
//                        + ":" + segments.get(
//                        segments.size() - 1).getTo());
//            }
            // List protein segments

            for (int i = 0; i < domains.size(); i++) {
                Domain domSeed = domains.get(i);
                Repeat repeat = new Repeat();
                repeat.getRepeats().add(getRepeatContent(domSeed));
                // right
                for (int j = i + 1; j < domains.size(); j++) {
                    Domain domCompare = domains.get(j);
                    if (isSameDomain(domSeed, domCompare, ss)) {
                        repeat.getRepeats().add(getRepeatContent(domCompare));
                    } else {
                        break;
                    }
                }
                // left
                for (int j = i - 1; j > 0; j--) {
                    Domain domCompare = domains.get(j);
                    if (isSameDomain(domSeed, domCompare, ss)) {
                        repeat.getRepeats().add(getRepeatContent(domCompare));
                    } else {
                        break;
                    }
                }
                if (repeat.getRepeats().size() >= 2) {
                    repeat.sortByPosition();
                    repeats.add(repeat);
                }
            }



            this.combineScore.setTmScore(maxTM);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());

        }
    }


    private RepeatContent getRepeatContent(Domain domain) {
        List<Segment> segments = domain.getSegments();
        return new RepeatContent(segments.get(0).getFrom(), segments.get(
                segments.size() - 1).getTo());
    }


    private List<Domain> removeOverlapDomain(List<Domain> domains) {

        List<Domain> newDomains = new ArrayList<Domain>();

        for (Domain domain : domains) {

            boolean isOverlap = false;
            for (Domain d : newDomains) {
                if (this.isOverlapDomain(domain, d)) {
                    isOverlap = true;
                    break;
                }

            }

            if (!isOverlap) {
                newDomains.add(domain);
            }
        }

        return newDomains;
    }

    public boolean isOverlapDomain(Domain domain1, Domain domain2) {

        List<Segment> segments1 = domain1.getSegments();
        List<Segment> segments2 = domain2.getSegments();

        // check with overlapping

        int start1 = segments1.get(0).getFrom();
        int end1 = segments1.get(
                segments1.size() - 1).getTo();

        int start2 = segments2.get(0).getFrom();
        int end2 = segments2.get(
                segments2.size() - 1).getTo();

        if ((start1 < start2 && start2 < end1) || (start1 < end2 && end2 < end1)) {
            return true;
        }
        return false;
    }

    public boolean isSameDomain(Domain domain1, Domain domain2, String strSS)
            throws StructureException {

        List<Segment> segments1 = domain1.getSegments();
        List<Segment> segments2 = domain2.getSegments();

        // check with overlapping

        int start1 = segments1.get(0).getFrom();
        int end1 = segments1.get(
                segments1.size() - 1).getTo();

        int start2 = segments2.get(0).getFrom();
        int end2 = segments2.get(
                segments2.size() - 1).getTo();

        if ((start1 < start2 && start2 < end1) || (start1 < end2 && end2 < end1)) {
            return false;
        }


        String pattern1 = VectorShape.getSSPattern(strSS
                .substring(start1, end1));
        String pattern2 = VectorShape.getSSPattern(strSS
                .substring(start2, end2));
        Atom[] atomSet1 = convertDomainToAtomSet(domain1);
        Atom[] atomSet2 = convertDomainToAtomSet(domain2);
//        boolean check = true;
//        if (atomSet1.length < 25) {
//            check = pattern1.equals(pattern2) && (pattern1.length() > 1 || pattern1.length() == 0);
//        } else if (atomSet1.length < 40) {
//            if (pattern1.length() > 1 && pattern1.length() <= 3)
//                check = pattern1.equals(pattern2);
//            else if (pattern1.length() == 1 && pattern1.equals(pattern2))
//                check = false;
//            else
//                check = true;
//        } else if (pattern1.length() == 1 && pattern1.equals(pattern2)) {
//            check = false;
//        }

        boolean check = (pattern1.equals(pattern2) && pattern1.length() == 2) || (pattern1.length() >= 2 && pattern2.length() >= 2);
        if(check) {
            SuperimposerOutput superOutput = superimposer.compareTwoStructures(atomSet1, atomSet2);
            maxTM = Math.max(maxTM, superOutput.getTmScore());
            return superOutput.isTRs();
        }

        return false ;

    }

    private Atom[] convertDomainToAtomSet(Domain domain) {
        List<Atom> listAtoms = new ArrayList<Atom>();
        List<Segment> segments = domain.getSegments();
        for (Segment s : segments) {
            // add
            Atom[] tmpatoms = Fragement.getFragementsofAtoms(atoms,
                    s.getFrom(), s.getTo());

            for (Atom a : tmpatoms) {
                listAtoms.add(a);
            }
        }
        return listAtoms.toArray(new Atom[listAtoms.size()]);
    }


}
