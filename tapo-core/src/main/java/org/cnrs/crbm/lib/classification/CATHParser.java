package org.cnrs.crbm.lib.classification;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathFactory;
import org.biojava.nbio.structure.cath.CathInstallation;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.util.List;

/**
 * Created by pdoviet on 10/5/2015.
 */
public class CATHParser {


    public static void main(String[] args) {
        String pdbID = "4ind";
        String pdbChain = "A";
        CathInstallation cath = new CathInstallation(Dir.CATH_LOCAL);
        cath.setCathVersion("4.1.0");

        try {
            List<CathDomain> domains = cath.getDomainsForPdb(pdbID);
            // show the structure in 3D
            for (CathDomain domain : domains) {
                System.out.println(domain.getCATH());
                System.out.println(domain.getRanges());
            }
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    CathInstallation cath;

    public CATHParser() {
        cath = new CathInstallation(Dir.CATH_LOCAL);
        cath.setCathVersion("4.1.0");

    }

    public String exportProperties(String fileData) {

        // read the pdbs from CATH v.4.1.0
        List<String> pdbs = DataIO.readLines(fileData);

        StringBuilder builder = new StringBuilder();
        // builder.append("pdbCode,pdbChain,scope,H,B,HH,HB,BH,BB\n");
//        ProgressBar bar = new ProgressBar();
//        int process = 1;
//        int sizeOfProcess = pdbs.size();

        for (String pdbID : pdbs) {
//            bar.update(process, sizeOfProcess);
//            process++;
            List<CathDomain> domains = cath.getDomainsForPdb(pdbID);
            try {
                // show the structure in 3D
                for (CathDomain domain : domains) {
//                    System.out.println(domain.getCATH());

                    for (String rangeCATH : domain.getRanges()) {
                        String pdbChainCATH = rangeCATH.split("_")[0];
                        RepeatFinder repeatFinder = new RepeatFinder(pdbID, pdbChainCATH);
                        //String[] rs = rangeCATH.split("_")[1].split("-");
                        Atom[] atoms = repeatFinder.getAtoms();

                        int startRef = 0;
                        int endRef = 0;
                        // convert to reference.
                        String region = rangeCATH.split("_")[1];
                        if (region.startsWith("-")) {
                            region = region.substring(1, region.length());
                            startRef = (-1) * Integer.parseInt(region.split("-")[0]);
                            endRef = Integer.parseInt(region.split("-")[1]);

                        } else {
                            startRef = Integer.parseInt(region.split("-")[0]);
                            endRef = Integer.parseInt(region.split("-")[1]);
                        }

                        startRef = PdbTools.getPosition(atoms, startRef);
                        endRef = PdbTools.getPosition(atoms, endRef);

                        String ssStr = repeatFinder.getStrSS().substring(startRef, endRef + 1);
                        String ssPattern = VectorShape.getSSPattern(ssStr);
                        double perH = ScopeParser.percentOfSecondaryStructureByType(ssStr, 'H');
                        double perB = ScopeParser.percentOfSecondaryStructureByType(ssStr, 'B');
                        double perHH = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "HH");
                        double perHB = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "HB");
                        double perBH = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "BH");
                        double perBB = ScopeParser.percentageArrangemanceOfStructure(ssPattern, "BB");
                        builder.append(pdbID + "_" + pdbChainCATH + "," + domain.getCATH() + "," + ssStr.length() + ",");
                        builder.append(NumberFormatUtils.format(perH) + ",");
                        builder.append(NumberFormatUtils.format(perB) + ",");
                        builder.append(NumberFormatUtils.format(perHH) + ",");
                        builder.append(NumberFormatUtils.format(perHB) + ",");
                        builder.append(NumberFormatUtils.format(perBH) + ",");
                        builder.append(NumberFormatUtils.format(perBB) + "\n");

//                        System.out.println(builder.toString());

                        //if(range.length())
                    }


                }
            } catch (Exception e) {
                // TODO Auto-generated catch block
                //e.printStackTrace();
            }

        }

        return builder.toString();


    }


}
