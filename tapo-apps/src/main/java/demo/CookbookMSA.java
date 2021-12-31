package demo;


import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.util.ConcurrencyTools;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;

public class CookbookMSA {

    public static void main(String[] args) throws Exception {



        List<ProteinSequence> lst = new ArrayList<ProteinSequence>();
        ProteinSequence seq1 = new ProteinSequence("MID");
        lst.add(seq1);
        ProteinSequence seq2 = new ProteinSequence("PTR");
        lst.add(seq2);
        ProteinSequence seq3 = new ProteinSequence("QFC");
        lst.add(seq3);
        ProteinSequence seq4 = new ProteinSequence("QFC");
        lst.add(seq4);
        ProteinSequence seq5 = new ProteinSequence("QFC");
        lst.add(seq5);
        ConcurrencyTools.setThreadPoolSingle();

        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        //System.out.printf("Clustalw:%n%s%n", profile);
        String msa = profile.toString();
        System.out.println(msa);
        //System.out.println(msa.split("\n")[0]);

        ConcurrencyTools.shutdown();

        //new CookbookMSA().test();


//        List<SiftsEntity> entities = SiftsMappingProvider.getSiftsMapping("1gix");
//        for (SiftsEntity e : entities){
//            System.out.println(e.getEntityId() + " " +e.getType());
//
//            for ( SiftsSegment seg: e.getSegments()) {
//                System.out.println(" Segment: " + seg.getSegId() + " " + seg.getStart() + " " + seg.getEnd()) ;
//
//                for ( SiftsResidue res: seg.getResidues() ) {
//                    System.out.println("  " + res);
//                }
//            }
//
//        }
//		String[] ids = new String[] { "Q21691", "Q21497", "O48771" };
//		try {
//			multipleSequenceAlignment(ids);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
    }

    private void test() {

        Manager manager = new Manager();
        manager.setName("Andrey");
        manager.setRole("superviser");
        Person person = manager;
        Manager manager2 = new Manager();
        manager2.setName("Phuong");
        manager2.setRole("student");


        System.out.println(person);
        System.out.println(manager2);
    }


    public abstract class Person {
        public String getName() {
            return name;
        }

        public void setName(String name) {
            this.name = name;
        }

        private String name;
    }

    public class Manager extends Person {
        public String getRole() {
            return role;
        }

        public void setRole(String role) {
            this.role = role;
        }

        private String role;

        public String toString() {
            return this.getName() + " " + this.getRole();

        }

    }

    private static void multipleSequenceAlignment(String[] ids)
            throws Exception {
        List<ProteinSequence> lst = new ArrayList<ProteinSequence>();
        for (String id : ids) {
            lst.add(getSequenceForId(id));
        }
        Profile<ProteinSequence, AminoAcidCompound> profile = Alignments
                .getMultipleSequenceAlignment(lst);
        System.out.printf("Clustalw:%n%s%n", profile);
        ConcurrencyTools.shutdown();
    }

    private static ProteinSequence getSequenceForId(String uniProtId)
            throws Exception {
        URL uniprotFasta = new URL(String.format(
                "http://www.uniprot.org/uniprot/%s.fasta", uniProtId));
        ProteinSequence seq = FastaReaderHelper.readFastaProteinSequence(
                uniprotFasta.openStream()).get(uniProtId);
        System.out.printf("id : %s %s%n%s%n", uniProtId, seq,
                seq.getOriginalHeader());
        return seq;
    }

}