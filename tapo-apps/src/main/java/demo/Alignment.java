package demo;

/**
 * Created by pdoviet on 7/1/2015.
 */

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.cnrs.crbm.lib.repeats.module.TMEvaluation;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

public class Alignment {

    public static void main(String[] args) throws CompoundNotFoundException {

//        AlignmentGui.getInstance();


//       / System.out.println("DNKNFYFRNGLPQIGVFK".substring(1, 3 + 1));

//        String[] copies = new String[]{"TAIGANGYKIIDNKNFYFRNGLPQIGVFK", "GPNGFEYFAPANTDANNIEGQAIRYQNRFLH", "LLGNIYYFGNNSKAVTGWQT", "INGNMYYFMPDTAMAAAGGLFE", "IDGVIYFFGVDGVKA"};
//        // Create multialignment.
//        ProteinSequence p1 = new ProteinSequence(copies[0]);
//        p1.setAccession(new AccessionID("seq1"));
//        ProteinSequence p2 = new ProteinSequence(copies[1]);
//        p2.setAccession(new AccessionID("seq2"));
//        ProteinSequence p3 = new ProteinSequence(copies[2]);
//        p3.setAccession(new AccessionID("seq3"));
//        ProteinSequence p4 = new ProteinSequence(copies[3]);
//        p4.setAccession(new AccessionID("seq4"));
//        ProteinSequence p5 = new ProteinSequence(copies[4]);
//        p5.setAccession(new AccessionID("seq5"));
//        ArrayList<ProteinSequence> lst = new ArrayList<ProteinSequence>();
//        lst.add(p1);
//        lst.add(p2);
//        lst.add(p3);
//        lst.add(p4);
//        lst.add(p5);
//
//
//        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
//        System.out.printf("Clustalw:%n%s%n", profile);
//        ConcurrencyTools.shutdown();
////
//        List<AlignedSequence<ProteinSequence, AminoAcidCompound>> alSeq = profile.getAlignedSequences();
//        System.out.println(profile.toString());
//        for (Sequence<AminoAcidCompound> seq : alSeq) {
//            ProteinSequence pSeq = new ProteinSequence(seq.getSequenceAsString());
//            pSeq.setAccession(seq.getAccession());
//            System.out.println(pSeq.getAccession() + "" + pSeq.getSequenceAsString());
//
//        }


//        Msa ms = new Msa(copies, "");
//        LinkedList<String> aligned = ms.buildAlignment();
//        for (String a : aligned) {
//            System.out.println(a);
//        }

        TMEvaluation tmEvaluation = new TMEvaluation();
        String[] copies = new String[]{"HB", "BH"};

        System.out.println(tmEvaluation.scoreSSPattern(copies));
        copies = new String[]{"HBB", "BHB"};
        System.out.println(tmEvaluation.scoreSSPattern(copies));


//        System.out.println("ABCDEFGHIJK".substring(1,4));
//        String[] ids = new String[] {"Q21691", "Q21691"};
//        try {
//            multipleSequenceAlignment(ids);
//        } catch (Exception e){
//            e.printStackTrace();
//        }
    }

    private static void multipleSequenceAlignment(String[] ids) throws Exception {
        List<ProteinSequence> lst = new ArrayList<ProteinSequence>();
        for (String id : ids) {
            lst.add(getSequenceForId(id));
        }
        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(lst);

        System.out.printf(":%n%s%n", profile);
        //ConcurrencyTools.shutdown();
    }

    private static ProteinSequence getSequenceForId(String uniProtId) throws Exception {
        URL uniprotFasta = new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uniProtId));
        ProteinSequence seq = FastaReaderHelper.readFastaProteinSequence(uniprotFasta.openStream()).get(uniProtId);
        System.out.printf("id : %s %s%n%s%n", uniProtId, seq, seq.getOriginalHeader());
        return seq;
    }

}