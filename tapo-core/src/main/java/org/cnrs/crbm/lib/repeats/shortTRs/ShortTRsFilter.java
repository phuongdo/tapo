package org.cnrs.crbm.lib.repeats.shortTRs;

import org.cnrs.crbm.lib.trsfinder.Repeat;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 8/30/2015.
 *
 * from version 1.1.0
 *
 */

public class ShortTRsFilter {


    public static void main(String[] args) {
        //test
    }

    public static boolean filter(Repeat repeat, String strAp) {

        boolean check = true;
        int L = repeat.getEnd() - repeat.getStart() + 1;
        int halfL = L / 2;

        if(L<0){
            System.out.println();
        }
        if (L < 40) {

            String regex = ".*[agfh]{" + halfL + ",}.*";
            if (strAp.matches(regex))
                check = false;

            regex = ".*[bv]{" + halfL + ",}.*";
            if (strAp.matches(regex))
                check = false;
        }
//
//        if (L < 40) {
//            for (int i = repeat.getStart(); i < repeat.getEnd() - halfL - 1; i++) {
//                String pattern = strAp.substring(i, i + halfL);
//
//
//            }
//
//
//        }




        return check;

    }


}
