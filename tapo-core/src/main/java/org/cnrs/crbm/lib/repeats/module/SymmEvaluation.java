package org.cnrs.crbm.lib.repeats.module;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.trsfinder.Repeat;

import java.util.List;

/**
 * Created by pdoviet on 5/26/2015.
 */
public class SymmEvaluation {

    MutilAlign mutilAign = new MutilAlign();
    Superimposer superimposer = new Superimposer();

    public Repeat findTheBestRepeatAmongCandidates(List<Repeat> candidates, Features features) throws StructureException {

        Repeat theBest = new Repeat();
        double bscore = 0.0;
        for(Repeat candidate :candidates ){
            double score = this.calsScore(candidate,features);
            if(bscore<score){
                bscore = score;
                theBest = candidate;
            }
        }
        return theBest;

    }

    public double calsScore(Repeat candidate, Features features) {

        double score = 0.0;
        Atom[] atoms = features.getAtoms();
        for (int i = 0; i < candidate.getRepeats().size() - 1; i++) {
            try {
                Atom[] query = Fragement.getFragementsofAtoms(atoms, candidate.getRepeats().get(i).getStart(), candidate.getRepeats().get(i).getEnd());
                Atom[] target = Fragement.getFragementsofAtoms(atoms, candidate.getRepeats().get(i + 1).getStart(), candidate.getRepeats().get(i + 1).getEnd());
                AFPChain afpChain = mutilAign.pairAlign(query, target);
                score = Math.max(afpChain.getTMScore(),score);
            } catch (Exception ex) {
                ex.printStackTrace();
            }

        }

        return score;

    }
}
