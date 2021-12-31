package org.cnrs.crbm.lib.trsfinder;

import org.cnrs.crbm.lib.repeats.module.HetatmModule;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

/**
 * Created by pdoviet on 5/26/2015.
 */
public class LigandBaseFinder extends Finder {
    static Logger logger = LoggerFactory.getLogger(LigandBaseFinder.class);

    public LigandBaseFinder(Features features) {
        super(features);
        this.name = "LigandBase";
    }

    @Override
    public void findRepeat(Features features) {

        try {
            HetatmModule hetatmModule = new HetatmModule();
            List<Repeat> lres = hetatmModule.findRepeats(features);
            this.repeats.addAll(lres);
            combineScore.setlScore(hetatmModule.getMaxHScore());

        } catch (Exception e) {
            e.printStackTrace();
            logger.error(e.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }


    }
}

