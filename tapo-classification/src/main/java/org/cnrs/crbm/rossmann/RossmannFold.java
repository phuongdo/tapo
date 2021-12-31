package org.cnrs.crbm.rossmann;


import org.biojava.nbio.structure.scop.ScopInstallation;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.repeats.CombineScore;
import org.cnrs.crbm.trclassification.ProteinSCOP;
import org.cnrs.crbm.trclassification.ScopeClassifierML;

import java.util.Map;

/**
 * Created by pdoviet on 7/12/2015.
 * Classify and analyse rossmann fold.
 */
public class RossmannFold {

    public static void main(String[] args) throws Exception {
        RossmannFold rossmannFold = new RossmannFold();
        rossmannFold.process();
    }

    private void extractDataToCSV() throws Exception {

        // read the output from TAPO analyses.(rossmann.TRs.o)
        ReadFasta fasta = new ReadFasta();
        Map<String, String> a = fasta.getFastaFileInSpecialCase("data/rossmann/rossmann.TRs.o");
        // load scope trclassification
        String cacheLocation = Dir.SCOPE_LOCAL;

        // download SCOP if required and load into memory
        ScopInstallation scop = new ScopInstallation(cacheLocation);

        for (Map.Entry<String, String> entry : a.entrySet()) {
//            System.out.println(entry.getKey());
            //System.out.println(">>>>>>" + entry.getValue());
            String pdb = entry.getKey().split("\\|")[0];
            CombineScore combineScore = new CombineScore(entry.getKey().split("\\|")[1]);
            System.out.println(pdb + "," + combineScore.toCSVFormat());


        }
    }

    private void process() throws Exception {

        // read the output from TAPO analyses.(rossmann.TRs.o)
        ReadFasta fasta = new ReadFasta();
        Map<String, String> fastaTRs = fasta.getFastaFileInSpecialCase("data/rossmann/rossmann.TRs.o");
        // load scope trclassification
        String cacheLocation = Dir.SCOPE_LOCAL;
        //String content = fastaTRs.get("1zp3_A|TM-Score=0.567;Psim-Score=1.000;CE-Score=4.521;V-Score=0.543;L-Score=0.227;CA-Score=0.000;S-Score=0.681");


        System.exit(1);
        // download SCOP if required and load into memory
        ScopInstallation scop = new ScopInstallation(cacheLocation);

        for (Map.Entry<String, String> entry : fastaTRs.entrySet()) {
//            System.out.println(entry.getKey());
            //System.out.println(">>>>>>" + entry.getValue());
            String pdb = entry.getKey().split("\\|")[0];
            String pdbCode = pdb.substring(0, 4);
            String pdbChain = pdb.substring(5, 6);
            ProteinSCOP scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, 0, 0);

            if (entry.getValue().contains("No-TRs")) {
                System.out.println(pdb + "\t" + scopeDomain + "\t" + 0);
            } else {
                System.out.println(pdb + "\t" + scopeDomain + "\t" + 1);
            }

        }
    }


}
