package org.cnrs.crbm.lib.largescale;


import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopInstallation;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.Row;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 6/3/2015.
 */
public class BuildNegativeSet {
    public static void main(String[] args) throws Exception {

        String cacheLocation = Dir.SCOPE_LOCAL;
        // download SCOP if required and load into memory
        ScopInstallation scop = new ScopInstallation(cacheLocation);
        BuildNegativeSet negativeSet = new BuildNegativeSet();


        LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
                .readFastaProteinSequence(new File("data/resultLargeScale/TAPO_NOTRs_07.cdhit"));
        for (Map.Entry<String, ProteinSequence> entry : a.entrySet()) {

            Row row = new Row();
            row.setLength(entry.getValue().getLength());
            row.setProtein(entry.getValue().getOriginalHeader().substring(0, 6));

            String strClass = negativeSet.classifyProteinByScope(row.getPdbCode(), row.getPdbChain(), scop);
            System.out.println(row.getProtein() + "\t" + strClass + "\t" + row.getLength());

        }

    }

    public String classifyProteinByScope(String pdbCode, String pdbChain, ScopInstallation scop) {


        String pdbId = pdbCode.toUpperCase();
        String strClass = "unk";
        // scop = (ScopInstallation)
        // ScopFactory.getSCOP(ScopFactory.VERSION_1_75);
        List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);
        for (ScopDomain domain : domains) {

            if (domain.getChains().contains(pdbChain)) {
                //System.out.println("XXXX:" + domain.getClassId());

                String description = scop.getScopDescriptionBySunid(domain.getClassId()).toString();

                if (description.contains("All alpha proteins")) {
                    strClass = "a";
                } else if (description.contains("All beta proteins")) {
                    strClass = "b";
                } else if (description.contains("Alpha and beta proteins (a/b)")) {
                    strClass = "a/b";
                } else if (description.contains("Alpha and beta proteins (a+b)")) {
                    strClass = "a+b";
                }

                //System.out.println(pdbCode + "_" + pdbChain + "\t" + strClass);
                //ScopNode node = scop.getScopNode(domain.getClassId());
                //System.out.println(scop.getScopDescriptionBySunid(domain.getClassId()));

//                while (node != null) {
//
//                    System.out.println("This node: sunid:" + node.getSunid());
//                    System.out.println(scop.getScopDescriptionBySunid(node.getSunid()));
//                    node = scop.getScopNode(node.getParentSunid());
//                }

            }
        }

        return strClass;


    }
}
