package org.cnrs.crbm.lib.repeats.shortTRs;

import org.cnrs.crbm.lib.sadb.ConformationMatcher;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;

import java.util.*;
import java.util.regex.Pattern;

/**
 * Created by pdoviet on 1/26/2015.
 * INPUT: conformational alphabets
 * OUTPUT: List of Regions containning pre-matched pattern
 */
public class ShortTRsDetector {
    // example comes from 3iox_A
    // public static final String EXAMPLE_TEST = "xaaaaaaaaaaaaaaaaaaafaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaafaaaaaaaaaaaaaaaaaaaaaaaaaaaaaagaaaaaaaaaaaaaaaaaaaaaaaggggbpibpabpbpvgvbxabxagpbvbbbbgxbpvpaaffaaffgbeghaaaagglbvbfggsphaghgpaggpbgbpbpbkxbegxplpbvgaxegppvbvvvvbvbplpbvbvbvblsagvbbllppvabvbvbvbspagsgsapgbvvvvvsapfhgkvvvvhpxxbpvplvvvvvpvbvbvbpaglbpvbsxlsvvsvblxpvbggbvbvvblbpebpvpvplpgsbbbllbvvsagbsxabplbllsabvgbpxgpplplfxbagpaxgghxvvvvpvbegbvbvvvvvvsfggabpaggppbvplbbggppbpxbvbvbvbhxpgpblppbpgbpppbppbppbpppppbbbbpbpppppppppbpppppbppbx";
    public static final String EXAMPLE_TEST = "xxabpvsaavbbvpfaahvabvbbvbbspaapaafaplvavpvvbaaffabepxabbabbvsaaffaxvgsbbssplvsvbvapxbxlxfxbbvvvapppvvvvafsvsvvbpplavbfaffaaffasaafaaafgaaeaebxbxbvpaaafaasbvsvspaaffaaffaafaafaaaaaafaaaaaafgaafaaafaaffaakbxxbbbpbabepvsaasbbbbfaafvbvvbvbbspaapgafapxxvsvvbfaffaaepvabvabbvpaafaaaalsbpvpplsvvvvgblkxxesbvbvvvsbbbgabbvvvvafvvvvbsaxpbvaaffafffgsgafaaafaaavbpbapblbvpaaafagppvvvspaaffaaffaafaafaaafaafaaafaafaaaaaaffaaffaafhpafaappppbabgaaaapaaaaaafaaafaaaaaaafaabaaaaax";


    public static void main(String[] args) {


        ShortTRsDetector detector = new ShortTRsDetector();
        List<LocationMatch> list = detector.getRepeatLocations(EXAMPLE_TEST);
//        System.out.println(EXAMPLE_TEST.length());
//
//        for (LocationMatch match : list) {
//
//            System.out.println(match.getStart() + "-" + match.getEnd() + ":" + match.getMeg());
//        }

        List<Repeat> repeats = detector.convertToRepeats(list);
        for (Repeat repeat : repeats) {
            System.out.println(repeat);
        }


    }


    public List<Repeat> convertToRepeats(List<LocationMatch> list) {

        List<Repeat> repeats = new ArrayList<Repeat>();
        for (LocationMatch match : list) {
            //System.out.println(match.getStart() + "-" + match.getEnd() + ":" + match.getMeg());
            int legTR = 7;
            Repeat repeat = new Repeat();
            if (match.getMeg().equals("isLongHelix")) {
                // 7 residues repeat
                legTR = 7;
            } else if (match.getMeg().equals("isClass2")) {
                // 3 residues repeat
                legTR = 3;
            } else if (match.getMeg().equals("isClass2OnlyP")) {
                // 3 residues repeat
                legTR = 3;
            } else if (match.getMeg().equals("isClass1")) {
                // 1 residues repeat
                legTR = 3;
            } else if (match.getMeg().equals("case4")) {
                // 1 residues repeat
                legTR = 3;
            }

            for (int i = match.getStart(); i < match.getEnd() - legTR; i = i + legTR) {
                repeat.getRepeats().add(new RepeatContent(i, i + legTR - 1));
            }

            if (repeat.getRepeats().size() >= 2) {
                repeats.add(repeat);
            }

        }

        return repeats;

    }


    public static Map<Character, Integer> getCharFreq(String s) {
        Map<Character, Integer> charFreq = new HashMap<Character, Integer>();
        if (s != null) {
            for (Character c : s.toCharArray()) {
                Integer count = charFreq.get(c);
                int newCount = (count == null ? 1 : count + 1);
                charFreq.put(c, newCount);
            }
        }
        return charFreq;
    }


    class ValueComparator implements Comparator<Character> {

        Map<Character, Integer> base;

        public ValueComparator(Map<Character, Integer> base) {
            this.base = base;
        }

        // Note: this comparator imposes orderings that are inconsistent with equals.
        public int compare(Character a, Character b) {
            if (base.get(a) >= base.get(b)) {
                return -1;
            } else {
                return 1;
            } // returning 0 would merge keys
        }
    }

    public List<LocationMatch> getRepeatLocations(String seqletters) {

        List<LocationMatch> locationMatches = new ArrayList<LocationMatch>();

        /**
         * < 40 residues with 80% of similar alphabets
         */
        if (seqletters.length() < 50) {
            String lqs = seqletters.replaceAll("k", "e");
            lqs = lqs.replaceAll("i", "l");
            lqs = lqs.replaceAll("v", "b");
            lqs = lqs.replaceAll("s", "p");
            lqs = lqs.replaceAll("c", "d");
            lqs = lqs.replaceAll("h", "g");
            lqs = lqs.replaceAll("f", "a");
            Map<Character, Integer> map = this.getCharFreq(lqs);
            //Sort keys by values.
            ValueComparator bvc = new ValueComparator(map);
            TreeMap<Character, Integer> sorted_map = new TreeMap<Character, Integer>(bvc);
            sorted_map.putAll(map);
            Integer count = sorted_map.firstEntry().getValue();
            //System.out.println(sorted_map.firstEntry().getKey());
            char topChar = sorted_map.firstEntry().getKey();
            if ((double) count / seqletters.length() >= 0.8) {

                //build regex
                String regex = ".*[" + topChar + "]{" + (count / 2) + ",}.*";
                Pattern pattern = Pattern.compile(regex);
                if (lqs.matches(regex)) {
                    locationMatches.add(new LocationMatch(0, seqletters.length() - 1, "case4"));
                }

            }

        } else {


            // searching for long helix
            int winsize = 40;
//            int minStart = 100000;
//            int maxEnd = -1;

            int[] labels = new int[seqletters.length()];


            for (int i = 0; i < seqletters.length() - winsize; i++) {
                if (ConformationMatcher.isLongHelix(seqletters.substring(i, i + winsize))) {
//                    minStart = Math.min(i, minStart);
//                    maxEnd = Math.max(maxEnd, i + winsize);
                    for (int j = i; j < i + winsize; j++) {
                        labels[j] = 1;
                    }

                }
            }


            List<Region> lstRegion = this.findRegion(labels);
            for (Region r : lstRegion) {
                locationMatches.add(new LocationMatch(r.getStart(), r.getEnd(), "isLongHelix"));
            }
//            System.out.println(seqletters.substring(minStart,maxEnd));


//            minStart = 100000;
//            maxEnd = -1;

            // is Class1
            winsize = 20;
            labels = new int[seqletters.length()];
            for (int i = 0; i < seqletters.length() - winsize; i++) {

                if (ConformationMatcher.isClass1(seqletters.substring(i, i + winsize))) {
                    for (int j = i; j < i + winsize; j++) {
                        labels[j] = 1;
                    }
                }
            }


            lstRegion = this.findRegion(labels);
            for (Region r : lstRegion) {
                locationMatches.add(new LocationMatch(r.getStart(), r.getEnd(), "isClass1"));
            }
            // is Class2

            labels = new int[seqletters.length()];
            for (int i = 0; i < seqletters.length() - winsize; i++) {

                if (ConformationMatcher.isClass2(seqletters.substring(i, i + winsize))) {
                    for (int j = i; j < i + winsize; j++) {
                        labels[j] = 1;
                    }
                }
            }

            lstRegion = this.findRegion(labels);
            for (Region r : lstRegion) {
                locationMatches.add(new LocationMatch(r.getStart(), r.getEnd(), "isClass2"));
            }


            // is Class2Only P
            labels = new int[seqletters.length()];
            for (int i = 0; i < seqletters.length() - winsize; i++) {

                if (ConformationMatcher.isClass2OnlyP(seqletters.substring(i, i + winsize))) {
                    for (int j = i; j < i + winsize; j++) {
                        labels[j] = 1;
                    }
                }
            }

            lstRegion = this.findRegion(labels);
            for (Region r : lstRegion) {
                locationMatches.add(new LocationMatch(r.getStart(), r.getEnd(), "isClass2OnlyP"));
            }


        }
        return locationMatches;

    }

    public List<Region> findRegion(int[] labels) {
        List<Region> lstRegions = new ArrayList<Region>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < labels.length; i++) {
            double ch = labels[i];
            if (ch == 1) {
                start = i;
                while (i < labels.length
                        && labels[i] == 1) {
                    i++;
                }
                i = end = i - 1;
                // save
                Region rc = new Region(start, end);
                lstRegions.add(rc);
            }


        }
        return lstRegions;
    }


}
