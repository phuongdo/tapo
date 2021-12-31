package org.cnrs.crbm.lib.classification;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.util.*;

/**
 * Created by pdoviet on 7/8/2015.
 */
public class ScopeParser {


    public Map<String, String> getScopeMap() {
        return scopeMap;
    }

    Map<String, String> scopeMap = new HashMap<String, String>();

    public static void main(String[] args) {

        ScopeParser scopeParser = new ScopeParser();

//        DataIO.writeToFile(scopeParser.analyse(), "data/resultLargeScale/scopeDB20.corr.fex");

//        StringBuffer buffer = new StringBuffer();
//
//        List<String> lines = DataIO.readLines("data/resultLargeScale/scopeDB20.fex");
//        for(String line: lines){
//
//            int count = StringUtils.countMatches(line, ",");
//            if(count==8){
//
//                buffer.append(line+"\n");
//            }
//        }


    }


    public static double percentOfSecondaryStructureByType(String ss, char type) {
        //double percentage = 0.0;
        int count = 0;
        for (int i = 0; i < ss.length(); i++) {
            if (type == ss.charAt(i)) {
                count++;
            }
        }

        return (double) count / ss.length();
    }

    public static double percentageArrangemanceOfStructure(String ssPattern, String type) {
        int count = 0;
        for (int i = 0; i < ssPattern.length() - 1; i++) {
            if (ssPattern.substring(i, i + 2).equals(type))
                count++;
        }

        if (ssPattern.length() - 1 <= 0)
            return 0.0;
        else
            return (double) count / (ssPattern.length() - 1);
    }

    public  int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {

            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;

        }
        return 0;
    }

    public void load() {
        // extract data
        // LOAD DB
        String cacheLocation = Dir.SCOPE_LOCAL;
        String scopeDir = cacheLocation + "/dir.cla.scope.2.05-stable.txt";
        //
        List<String> lines = DataIO.readLines(scopeDir);
        for (String line : lines) {
            String[] cols = line.split("\t");
            String pdbCode = cols[1];
            String domainScope = cols[2] + "\t" + cols[3];
            scopeMap.put(pdbCode, domainScope);
        }

    }

    public String analyse() {

//        List<String> lines = DataIO.readLines("data/resultLargeScale/scopeDB20.err.fex");
//
//        Set<String> set = new HashSet<String>();
//        for (String line : lines) {
//
//            set.add(line.substring(0, 4));
//
//        }

        StringBuilder builder = new StringBuilder();
        // builder.append("pdbCode,pdbChain,scope,H,B,HH,HB,BH,BB\n");
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = scopeMap.size();
        for (Map.Entry<String, String> entry : scopeMap.entrySet()) {
            bar.update(process, sizeOfProcess);
            process++;

            try {
                String pdbCode = entry.getKey();
                String domainScope = entry.getValue();
//                if(!set.contains(pdbCode))
//                    continue;
                // System.out.println(pdbCode + ":" + domainScope);
                String[] details = domainScope.split("\t");

                List<String> lstDomain = new ArrayList<String>();
                if (details[0].contains(",")) {

                    for (String str : details[0].split(","))
                        lstDomain.add(str);

                } else {
                    lstDomain.add(details[0]);
                }


                for (String strDomain : lstDomain) {
                    String pdbChain = strDomain.substring(0, 1);
                    String scopeClass = details[1].substring(0, 1);
                    RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                    Atom[] atoms = repeatFinder.getAtoms();
                    int start = 0;
                    int end = atoms.length - 1;
                    if (strDomain.length() > 2) {
                        String region = strDomain.split(":")[1];
                        start = Integer.parseInt(region.split("-")[0]);
                        end = Integer.parseInt(region.split("-")[1]);
                        //convert to relative values
                        start = this.getPosition(atoms, start);
                        end = this.getPosition(atoms, end);
                    }

                    String ssStr = repeatFinder.getStrSS().substring(start, end + 1);
                    String ssPattern = VectorShape.getSSPattern(ssStr);
                    double perH = this.percentOfSecondaryStructureByType(ssStr, 'H');
                    double perB = this.percentOfSecondaryStructureByType(ssStr, 'B');
                    double perHH = this.percentageArrangemanceOfStructure(ssPattern, "HH");
                    double perHB = this.percentageArrangemanceOfStructure(ssPattern, "HB");
                    double perBH = this.percentageArrangemanceOfStructure(ssPattern, "BH");
                    double perBB = this.percentageArrangemanceOfStructure(ssPattern, "BB");
                    builder.append(pdbCode + "," + strDomain + "," + scopeClass + ",");
                    builder.append(NumberFormatUtils.format(perH) + ",");
                    builder.append(NumberFormatUtils.format(perB) + ",");
                    builder.append(NumberFormatUtils.format(perHH) + ",");
                    builder.append(NumberFormatUtils.format(perHB) + ",");
                    builder.append(NumberFormatUtils.format(perBH) + ",");
                    builder.append(NumberFormatUtils.format(perBB) + "\n");
                }
            } catch (Exception ex) {
//                System.out.println(entry.getKey() + "" + entry.getValue());
//                ex.printStackTrace();
            }


        }

        return builder.toString();


    }

    public ScopeParser(String scopeDir) {


        List<String> lines = DataIO.readLines(scopeDir);
        for (String line : lines) {
            String[] cols = line.split("\t");
            String pdbCode = cols[1];
            String domainScope = cols[2] + "\t" + cols[3];
            scopeMap.put(pdbCode, domainScope);
        }


    }

    public ScopeParser() {

        this.load();
        // do something

    }
}
