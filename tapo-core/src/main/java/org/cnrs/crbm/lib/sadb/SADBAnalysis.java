package org.cnrs.crbm.lib.sadb;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.cnrs.crbm.lib.conf.Dir;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

public class SADBAnalysis {

    private String MEG_LONGHELIX = "CL4";
    private String MEG_CALSS2_COLLAGEN = "C2C";
    private String MEG_CALSS2_ONLY_P = "C2P";
    private String MEG_CALSS1 = "CL1";

    public static void main(String[] args) throws Exception {
        new SADBAnalysis().run();

        // String seq =
        // "xgaaaaaaaaaaaaagaaaaaaaaaaaaaaafaaaaaaaaaaaaafaaaaaaaaax"
        // .toLowerCase();
        // System.out.println(new SADBAnalysis().annotate(seq));
    }

    public void run() throws Exception {
        // Try with the FastaReaderHelper
        LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
                .readFastaProteinSequence(new File(Dir.SADB_LOCAL + "/sadb"));

        System.out.println(a.size());
        for (Entry<String, ProteinSequence> entry : a.entrySet()) {
            // System.out.println(entry.getValue().getOriginalHeader() + "="
            // + entry.getValue().getSequenceAsString());

            System.out.println(entry.getValue().getSequenceAsString());

            String seq = entry.getValue().getSequenceAsString();
            String pdb = entry.getKey().substring(0, 6);
            // System.out.println(pdb);
//			String annotation = this.annotate(seq);
//			if (annotation.length() > 0) {
//				System.out.println(pdb + ":" + annotation);
//			}

        }
    }

    public String annotate(String seq) {

        seq = process(seq).toLowerCase();
        // System.out.println(seq);
        StringBuffer buffer = new StringBuffer();
        if (ConformationMatcher.isLongHelix(seq)) {
            buffer.append(MEG_LONGHELIX + ";");
        }
        if (ConformationMatcher.isClass1(seq)) {
            buffer.append(MEG_CALSS1 + ";");
        }
        if (ConformationMatcher.isClass2(seq)) {
            buffer.append(MEG_CALSS2_COLLAGEN + ";");
        }
        if (ConformationMatcher.isClass2OnlyP(seq)) {
            buffer.append(MEG_CALSS2_ONLY_P + ";");
        }

        return buffer.toString();

    }

    private String process(String seqs) {

        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < seqs.length(); i++) {
            char c = seqs.charAt(i);
            if (c == 'k') {
                c = 'e';
            } else if (c == 'i') {
                c = 'l';
            } else if (c == 'v') {
                c = 'b';
            } else if (c == 's') {
                c = 'p';
            } else if (c == 'h') {
                c = 'g';
            } else if (c == 'f') {
                c = 'a';
            }

            buffer.append(c);

        }

        return buffer.toString();
    }
}
