package org.cnrs.crbm.nextgen;

import org.cnrs.crbm.lib.utils.FastaUtils;

import java.io.IOException;

/**
 * Created by pdoviet on 10/21/2015.
 */
public class ExportFastaFile {
    public static void main(String[] args) throws IOException {
        //FastaUtils.exportFastaFileFromListPDBRegion("data/classification/pdbNewTRs.txt", "data/classification/pdbNewTRs.fasta");
//        FastaUtils.exportFastaFileFromListPDB("data/nrdb/nrPDB26June2015_40.txt","data/nrdb/nrPDB26June2015_40.fasta");
        FastaUtils.exportFastaFileFromListPDB("data/nrdb/abc.txt","data/nrdb/abc.fasta");


    }
}
