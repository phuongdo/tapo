package org.cnrs.crbm.lib.largescale;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.utils.FastaUtils;

import java.io.File;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 6/2/2015.
 */
public class ExportFastaFile {
    public static void main(String[] args) throws Exception {


        FastaUtils.exportFastaFileFromListPDBRegion("data/phylo/tree.in", "data/phylo/tree.fasta");

    }
}
