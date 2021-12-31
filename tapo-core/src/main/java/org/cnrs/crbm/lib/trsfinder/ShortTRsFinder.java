package org.cnrs.crbm.lib.trsfinder;

import org.cnrs.crbm.lib.repeats.shortTRs.LocationMatch;
import org.cnrs.crbm.lib.repeats.shortTRs.ShortTRsDetector;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class ShortTRsFinder extends Finder {
    static Logger logger = LoggerFactory.getLogger(ShortTRsFinder.class);

    public ShortTRsFinder(Features features) {
        super(features);
        this.name = "ShortTRs";
    }

    @Override
    public void findRepeat(Features features) {
        try {
            ShortTRsDetector detector = new ShortTRsDetector();
            List<LocationMatch> listFromDetector = detector.getRepeatLocations(features.getStrSeqAp());
            List<Repeat> shortRepeats = detector.convertToRepeats(listFromDetector);
            this.repeats.addAll(shortRepeats);
            if (listFromDetector.size() > 0)
                this.combineScore.setCaScore(1);

        } catch (Exception ex) {
            logger.error(ex.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }

        //this.repeats = repeats_copy;

    }
}
