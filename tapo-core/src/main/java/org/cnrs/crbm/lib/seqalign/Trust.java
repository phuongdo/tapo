package org.cnrs.crbm.lib.seqalign;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

import nl.vu.cs.align.Assert;
import nl.vu.cs.align.FakeProgressBar;
import nl.vu.cs.align.ProgressBar;
import nl.vu.cs.align.TextProgressBar;
import nl.vu.cs.align.algorithm.Align;
import nl.vu.cs.align.algorithm.AlignAffine;
import nl.vu.cs.align.algorithm.AlignDDPMin;
import nl.vu.cs.align.algorithm.ComparableSequence;
import nl.vu.cs.align.algorithm.CompleteAlignSelfAffine;
import nl.vu.cs.align.algorithm.LocalAlignAffineGotoh;
import nl.vu.cs.align.algorithm.LocalAlignAffineNice;
import nl.vu.cs.align.algorithm.LocalInit;
import nl.vu.cs.align.algorithm.Profile;
import nl.vu.cs.align.algorithm.SelfLocalFill;
import nl.vu.cs.align.algorithm.Sequence;
import nl.vu.cs.align.algorithm.Trace;
import nl.vu.cs.align.matrix.AlignData;
import nl.vu.cs.align.matrix.AlignDataAffine;
import nl.vu.cs.align.matrix.AlignDataAffineGotoh;
import nl.vu.cs.align.matrix.Matrix;
import nl.vu.cs.align.proteinlibs.FastaDb;
import nl.vu.cs.align.proteinlibs.LengthComparator;
import nl.vu.cs.align.proteinlibs.ProteinLib;
import nl.vu.cs.align.repeats.EstimateEVDParams;
import nl.vu.cs.align.substtable.SubstitutionTable;

public class Trust {

	String blosumDir;
	String pattern;
	String name;
	int gapOpen = -1;
	int gapExtend = -1;
	private Collection repeats_output;
	String aACID = "";
	public List<RepeatSeq> getRepeats() {
		return repeatSqs;
	}

	List<RepeatSeq> repeatSqs = new ArrayList<RepeatSeq>();

	public static void main(String[] args) {
		String x1BIK = "xagpppbaegpgxpbpmbpbgaglbsmbpbhxblxbgxmaboangahpbgghggpppbaegpxxpmmvmoglbpmmshxblmhxwgboafoag";

		String x1a50_new = "OOFFHTNPKNXVTQMARTDRPVRVPKEVWLPFEEXKKKISRTNEHEHKKCDPXPFHEFNKQPFKNMGPGQPKHHHFKKISRNHHEHZRHVRTTNHHEKKVSWQTQPFKTPKEVSQVQXMPHXPHEEFHKLSQTFEENKNXTMWQHHKHEHFFFEGLSRTHKHMNVRMPEKTZVTFFHHEEFNKISRTHEHHFFEEVVQPDGQNZTHFHEFEMSVQVPEFKISRTKFNGSHKILQTKHVPHFKDQSRKFMRKHEFKKISRNHEFXNHKKKACPQPHKNHEHKETLLTHFEFEXTHEXPHKFQOOOOGSKTHEEEEVWRPHEFHKHNDDSTQTPEEEEEVPGLQTKHFHKEMPGTHHMQNKHEKNMBDALLSRTHEFHXNNFKOOOOKNTKNHKKFFFOO";
		// String strR2 =
		// "AESTLGAAAAQSGRYFGTAIASGKLGDSAYTTIASREFNMVTAENEMKIDATEPQRGQFNFSAGDRVYNWAVQNGKQVRGHTLAWHSQQPGWMQSLSGSTLRQAMIDHINGVMGHYKGKIAQWDVVNEAFSDDGSGGRRDSNLQRTGNDWIEVAFRTARAADPAAKLCYNDYNIENWTWAKTQGVYNMVRDFKQRGVPIDCVGFQSHFNSGSPYNSNFRTTLQNFAALGVDVAITELDIQGASSSTYAAVTNDCLAVSRCLGITVWGVRDTDSWRSGDTPLLFNGDGSKKAAYTAVLNALNGGSSTPPPSGGGQIKGVGSGRCLDVPNASTTDGTQVQLYDCHSATNQQWTYTDAGELRVYGDKCLDAAGTGNGTKVQIYSCWGGDNQKWRLNSDGSIVGVQSGLCLDAVGGGTANGTLIQLYSCSNGSNQRWTRT";
		int gapo = -8;
		int gape = -2;

		Trust trust = new Trust("BLOSUM62", "ABCDEFGHIKLMNPQRSTVWXYZO",
				"x1a50_new", x1a50_new.toUpperCase(), gapo, gape);
		Trust trust2 = new Trust("BLOSUM_COMBINED", "EKLIBVPSDCGHAFXYOWMN",
				"x1BIK", x1BIK.toUpperCase(), gapo, gape);

		try {
			System.out.println(trust.getRepeatsOutput().size());
			System.out.println(trust2.getRepeatsOutput().size());
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public Trust(String blosumDir, String aACID, String name, String pattern,
			int gapo, int gape) {

		this.blosumDir = blosumDir;
		this.pattern = pattern;
		this.name = name;
		this.gapExtend = gape;
		this.gapOpen = gapo;
		this.aACID = aACID;

		this.run();
	}

	public void run() {
		boolean PARAM_HTML = false;
		String[] arg = new String[0];

		SubstitutionTable blosum;
		ProteinLib library;
		ProteinLib librarySeg;
		// maximum protein length
		int gapo = this.gapOpen;
		int gapx = this.gapExtend;

		try {
			String libName = null;

			blosum = new SubstitutionTable(this.blosumDir, aACID);
			library = new FastaDb(this.name, this.pattern);
			librarySeg = library;

		} catch (IOException ioe) {
			throw new RuntimeException(ioe);
		}

		// cpu issues
		int procNr = 0;
		int procTotal = 1;

		// maximum protein length
		String outputDir = "pictures";
		int maxLen = Integer.MAX_VALUE;

		showAllTransitivity(library, librarySeg, maxLen, blosum, gapo, gapx,
				outputDir, procNr, procTotal, PARAM_HTML);

	}

	private static float prev_fitness;

	// constants
	public static final boolean VISUAL_PICT_INVERT = true;
	public static final boolean VISUAL_PICT_NORMAL = false;
	public static final boolean VISUAL_PROJECT_SW = false;
	public static final boolean VISUAL_PROJECT_SE = false;
	public static final boolean VISUAL_PROJECT_E = false;
	public static final boolean VISUAL_PROJECT_TRANS_ADD = false;
	public static final boolean VISUAL_PROJECT_NO_TRANS = false;
	public static final boolean VISUAL_PROJECT_CROSSSECTION = false;

	// logging constants

	// debug constants
	public static final boolean DEBUG_MUTUAL_ALIGN = false;
	public static final boolean DEBUG_VISUAL_CUTOFF = false; // the parameters
																// of the
																// alignments
																// found
	public static final boolean DEBUG_VISUAL_SIMPLE = false; // simple text
																// information
	public static final boolean DEBUG_VISUAL_EVERY_TRACE = false;
	public static final boolean DEBUG_VISUAL_SW = false; // "s-w" file in case
															// no trace found
	public static final boolean DEBUG_VISUAL_SUMMARY = false; // "profile wrapa",
																// "vector",
																// "trans add",
																// "result",
																// "s-w",
																// "SW proj"
	public static final boolean DEBUG_VISUAL_PARTIAL_TRACES = false;
	public static final boolean DEBUG_VISUAL_ESTIMATION = false; // parameters
																	// during
																	// statistical
																	// estimation:
																	// a, b,
																	// evdK,
																	// evdLambda
	public static final boolean DEBUG_VISUAL_STATS = false; // print all random
															// alignments
															// extracted
	public static final boolean DEBUG_REPORT_TIMES = false;
	public static final boolean DEBUG_PROFILE_TIMEING = false;
	public static final boolean DEBUG_PROFILE_SOURCE = false; // "profilesource"
																// file
	public static final boolean DEBUG_STATS_INSIGN = false; // show
															// multiplication of
															// insignificant
															// traces
	public static final boolean DEBUG_WE_UPDATE = false;
	public static final boolean DEBUG_WE_UPDATE_EVERY_TRACE = false; // shows
																		// the
																		// matrix
	public static final boolean DEBUG_WE_UPDATE_STATS = false;
	public static final boolean DEBUG_REPEAT_LEN_RETRIVAL = false;
	public static final boolean DEBUG_EXIT_AFTER_ESTIMATION = false; // exits
																		// after
																		// printing
																		// EVD
																		// params
	public static final boolean DEBUG_VECTOR_DECORATE = true; // thoroughly
																// decoration of
																// vector files
	public static final boolean DEBUG_PRINT_TRACES = false; // print
															// statistically
															// significant
															// traces and trans
															// add

	public static final boolean PARAM_SENSITIVE_TRACE = true;
	public static final int PARAM_SENSITIVE_TRACE_SMOOTH = 10;
	public static final int PARAM_SMOOTH_PROJECTION = 1;
	public static final int PARAM_WIDE_TRACE = 4; // inversion of the size of
													// trace from distance of
													// which it's treated as
													// noise
	public static final float PARAM_MIN_REPEAT = 0.5f; // minimal size of repeat
	/** filters traces around diagonal */
	public static final boolean PARAM_DIAG_FILTER = false; // TODO: noise is
															// introduced when
															// the filter
															// disabled unless
															// diagonal trace is
															// implemented
															// properly
	public static final int PARAM_MAX_ESTIM_SEQ_SIZE = 500; // for stat sign the
															// sequence size is
															// to be limited
	public static final int PARAM_MIN_TRACES_NUM = 3; // minimal number of
														// traces to extract
														// (possibly
														// insignificant, to be
														// validated later)
	public static boolean PARAM_FORCE = false;

	// assign gap penalties
	// private static final int gapo = -8;
	// private static final int gapx = -2;

	// TODO because X regions can be reported as repeats
	public Collection findRepeats(String seq, String seqOrig,
			SubstitutionTable subst, int gapo, int gapx, float pvalue,
			String proteinName, String outputDir, int typeNum, boolean html) {
		Assert.assertTrue(seq.length() == seqOrig.length());
		// This actually happens, sometimes you can find X in the sequence...
		// Assert.assertTrue(seqOrig.indexOf('x') < 0 && seqOrig.indexOf('X') <
		// 0);

		int maxIteration = 100;
		ProgressBar pb = null;
		if (DEBUG_REPORT_TIMES)
			pb = new TextProgressBar();
		else
			pb = new FakeProgressBar();

		pb.start("Proteins analyzed", 1000);

		pb.advance("Create AlignData");
		AlignAffine align = new LocalAlignAffineNice();

		// best alignments will be stored here
		Vector bestTraces = new Vector();
		Trace strongestTrace = null;
		double evdLambda = Double.NaN;
		double evdK = Double.NaN;

		float cutoff = -1f;
		// if (PARAM_PURE_CUTOFF) {
		// pb.info("Estimating cutoff score");
		// cutoff = 10.0f;
		// //System.err.println("Cutoff set to "+cutoff+" sequence independent!");
		// //cutoff = estimateCutoff(seq, subst);
		//
		// if (DEB UG_VISUAL_CUTOFF) {
		// System.out.println("cutoff "+cutoff);
		// }
		// }

		pb.advance("Starting extracting traces");
		long timeStart = System.currentTimeMillis();

		AlignDataAffine data = EstimateEVDParams.extractSignificantAlignments(
				new Sequence(seq), new Sequence(seq),
				true /* selfSimilarity */, align, subst, gapo, gapx, pvalue,
				maxIteration, bestTraces, false, /*
												 * no special treatment for
												 * repeats
												 */
				null, proteinName, typeNum, outputDir, pb);

		if (DEBUG_PRINT_TRACES) {
			// System.err.print("Vector traces");
			printTraces(seq, bestTraces);
		}

		strongestTrace = (Trace) (bestTraces.size() > 0 ? bestTraces
				.firstElement() : null);

		long timeEnd = System.currentTimeMillis();

		if (DEBUG_PROFILE_TIMEING)
			// System.err.println("Repeats extracted in " + (timeEnd -
			// timeStart));
			pb.advance("Repeats extracted in " + (timeEnd - timeStart));

		if (DEBUG_VISUAL_SUMMARY && bestTraces.size() == 0 && data != null) {
			// mirror
			// Matrix.mirrorTopRight(transM);
			// save the traces itself
			float m[][] = Matrix.copy(data.getMatrix());
			for (int y = 0; y < m.length; y++) {
				for (int x = 0; x < m[y].length; x++) {
					if (m[y][x] < 0)
						m[y][x] = 0;
				}
			}
			Matrix.saveImageGray(m, VISUAL_PICT_INVERT,
					createImageFileName(outputDir, proteinName, "s-w"));
		}

		// ------------------------------------------------------
		// here we assume that nicest repeat traces are extracted
		// ------------------------------------------------------

		// calculate the maximal tandem repeat size
		int maxTandemRepeatSize = strongestTrace != null ? strongestTrace
				.coveredMin() : 0;
		// calculateMaxTandemRepeatSize(bestTraces);

		/*
		 * if (VISUAL_BINARY_TRACE) { bestTraces = convertToBinary(bestTraces);
		 * }
		 */

		if (PARAM_SENSITIVE_TRACE) {
			// Assert.assertTrue(!VISUAL_BINARY_TRACE);
			bestTraces = convertToSensitive(bestTraces, data);
		} else {
			// System.err.println("Not converting the trace to sensitve");
		}

		// calculate the transitivity on "vector" traces
		pb.advance("Calculate transitive traces");
		Collection transitiveTraces = Trace.multiplyTraces(bestTraces);
		Trace.sort((List) transitiveTraces);

		if (DEBUG_PRINT_TRACES) {
			// System.err.print("Transitive traces");
			printTraces(seq, transitiveTraces);
		}

		// create the matrix to paint on
		pb.advance("Paint traces");

		float[][] transM = Matrix.createEmpty(data.getMatrix());
		// paint onto the matrix
		Trace.paintMax(transM, bestTraces);

		if (DEBUG_VISUAL_SUMMARY) {
			if (DEBUG_VECTOR_DECORATE) {
				saveImageDecorated(transM, proteinName, "vector " + typeNum,
						outputDir);
			} else {
				Matrix.saveImageGray(
						transM,
						VISUAL_PICT_INVERT,
						createImageFileName(outputDir, proteinName, "vector "
								+ typeNum));
			}

			if (VISUAL_PROJECT_NO_TRANS) {
				createHistograms(transM, proteinName, "no trans", outputDir);
			}
		}

		// TRANSITIVITY -- apply
		Trace.paintAdditively(transM, transitiveTraces);

		Matrix.clearDiagonal(transM);

		if (DEBUG_VISUAL_SUMMARY) {
			if (VISUAL_PROJECT_TRANS_ADD) {
				createHistograms(transM, proteinName, "trans add " + typeNum,
						outputDir);
			}
			// mirror
			// Matrix.mirrorTopRight(transM);
			// save the traces itself
			if (DEBUG_VECTOR_DECORATE) {
				saveImageDecorated(transM, proteinName, "vectortrans add "
						+ typeNum, outputDir);
			} else {
				Matrix.saveImageGray(
						transM,
						VISUAL_PICT_INVERT,
						createImageFileName(outputDir, proteinName,
								"vectortrans add"));
			}
		}

		// restrict to the area defined by the strongest trace
		int startMatrix = 1;
		int endMatrix = seq.length();
		int distanceFromDiagonal = 0;
		if (strongestTrace != null) {
			// TODO: In the rectangle shouldn't be defined necessarily around
			// diagonal
			// it should happen anywhere in the protein
			startMatrix = Math.min(strongestTrace.start[0],
					strongestTrace.end[0]);
			endMatrix = Math.max(
					strongestTrace.start[strongestTrace.length() - 1],
					strongestTrace.end[strongestTrace.length() - 1]);
			// TODO Heuristics: stack all traces one onto another and take the
			// boundaries of the strongest 90%. Currently: the strongest trace
			// is taken,
			// and it happens, that it's way too much
			distanceFromDiagonal = strongestTrace
					.medianDistanceToDiagonal(data);
			Assert.assertTrue(strongestTrace.start[0] <= strongestTrace.start[strongestTrace
					.length() - 1]);
			Assert.assertTrue(strongestTrace.end[0] <= strongestTrace.end[strongestTrace
					.length() - 1]);
		}
		Assert.assertTrue(endMatrix >= startMatrix);
		Assert.assertTrue(startMatrix > 0);
		Assert.assertTrue(endMatrix > 0);

		// int tandemRepeatSize = calculateTandemRepeatSize(transM,
		// maxTandemRepeatSize, startMatrix, endMatrix, outputDir, proteinName,
		// typeNum);
		Collection tandemRepeatSizes = calculateTandemRepeatSize(transM,
				maxTandemRepeatSize, startMatrix, endMatrix, outputDir,
				proteinName, typeNum);

		float[][] symTransM = transM;
		Matrix.mirrorTopRight(symTransM);
		Matrix.clearDiagonal(symTransM);
		transM = null; // sanity

		// best results out of different tandemRepeatSizes
		float bestRepeatsScore = -1; // anything is better than that...
		Collection bestRepeats = new Vector(); // return value has to be != null
		Profile bestProfile = null;
		double[] bestEvdParams = null;
		int bestStartWrapping = -1;

		// for every possible tandem repeat size, calculate the one leading to
		// the highest score
		Sequence sequence = new Sequence(seq);
		for (Iterator iter = tandemRepeatSizes.iterator(); iter.hasNext();) {
			int tandemRepeatSize = ((Integer) iter.next()).intValue();

			// OK, we have the offset of the trace: maxTraceOffset

			// do we have the size of tandem repeat or just the distance from
			// diagonal?

			// look for the highest scoring protein region of the length
			// maxTraceOffset
			pb.advance("Find the strap");
			int strapStart = maxScoringStrap(symTransM, startMatrix, endMatrix,
					tandemRepeatSize, distanceFromDiagonal);
			Assert.assertTrue(strapStart >= startMatrix);
			Assert.assertTrue(strapStart + tandemRepeatSize - 1 <= endMatrix);

			// correct the profile size

			float debugM[][] = DEBUG_PROFILE_SOURCE ? Matrix
					.createEmpty(symTransM) : null;

			pb.advance("Create the profile");
			Profile profile = null;
			float maxCol = Matrix.getMax(symTransM);

			if (tandemRepeatSize > 0) {
				// construct a profile of the length
				// <code>tandemRepeatSize</code>
				// TODO correct profile to include the diagonal!, of the
				// strength of the avg match!
				profile = new Profile(symTransM, new Sequence(seq),
						startMatrix, endMatrix - startMatrix + 1, strapStart,
						tandemRepeatSize, distanceFromDiagonal, subst, debugM);

				if (VISUAL_PROJECT_CROSSSECTION) {
					float m[] = Matrix.copy(symTransM[strapStart
							+ tandemRepeatSize / 2]);
					float maxV = Matrix.getMax(m);
					m[0] = maxV;
					float[][] histM = Matrix.paintHistogram(m, null, 100);
					Matrix.saveImageGray(
							histM,
							VISUAL_PICT_INVERT,
							createImageFileName(outputDir, proteinName,
									"clean proj"));
					m = Matrix.copy(symTransM[strapStart
							+ (int) (tandemRepeatSize * 1.5)]);
					m[0] = maxV;
					histM = Matrix.paintHistogram(m, null, 100);
					Matrix.saveImageGray(
							histM,
							VISUAL_PICT_INVERT,
							createImageFileName(outputDir, proteinName,
									"dirty proj"));
				}

				// draw lines on the matrix
				// Matrix.drawHLine(symTransM, startMatrix, strapStart,
				// endMatrix, strapStart, maxCol);
				// Matrix.drawHLine(symTransM, startMatrix,
				// strapStart+tandemRepeatSize-1, endMatrix,
				// strapStart+tandemRepeatSize-1, maxCol);
				// Matrix.drawVLine(symTransM, strapStart, startMatrix,
				// strapStart, endMatrix, maxCol);
				// Matrix.drawVLine(symTransM, strapStart+tandemRepeatSize-1,
				// startMatrix, strapStart+tandemRepeatSize-1, endMatrix,
				// maxCol);

			}

			if (DEBUG_VISUAL_SUMMARY) {
				if (DEBUG_VECTOR_DECORATE) {
					saveImageDecorated(symTransM, proteinName, "result "
							+ typeNum, outputDir);
				} else {
					Matrix.saveImageGray(
							symTransM,
							VISUAL_PICT_INVERT,
							createImageFileName(outputDir, proteinName,
									"result " + typeNum));
				}
				if (debugM != null) {
					Matrix.saveImageGray(
							debugM,
							VISUAL_PICT_INVERT,
							createImageFileName(outputDir, proteinName,
									"profilesource " + typeNum));
				}
			}

			// symTransM = null; // needed if the matrix has been modified

			pb.advance("Align against the profile");
			// bestTraces = null;

			// how the profile reflects on sequence
			Collection repeats = new Vector();
			data = null;
			if (profile != null && bestTraces.size() > 0) {

				boolean wrapa = true;
				double evdParams[] = new double[2];

				// TODO extract repeats based on non-filtered sequence
				AlignData dataSign = EstimateEVDParams
						.extractSignificantAlignments(profile, sequence,
								false /* selfSimilarity */,
								new LocalAlignAffineNice(), subst, gapo, gapx,
								pvalue, maxIteration, repeats, wrapa, /*
																	 * special
																	 * treatment
																	 * for the
																	 * case of
																	 * repeats
																	 * extraction
																	 */
								evdParams, proteinName, typeNum, outputDir, pb);
				// we have wrap-around repeats extracted

				if (repeats.size() > 0) {
					// convert wrap-around repeat into normal repeat
					int startWrapping = calcIntRepeatsFromWrapa(repeats);
					repeats = Trace.wrapaToNormal(startWrapping, repeats);

					// remove too short repeats
					// //System.err.println("remove short repeats removal (weak is enough)");
					// nonWrapa = filterTooShortRepeats(nonWrapa,
					// profile.length());

					// remove too weak repeats
					// System.err
					// .println("running into danger of removing everything!");
					// TODO remove repeats based on non-filtered sequence
					repeats = filterTooWeakRepeats(repeats, profile, sequence,
							subst, gapo, gapx, pvalue, evdParams[0],
							evdParams[1], proteinName);

					if (repeats.size() == 1) {
						// System.err
						// .println("!!! PROBLEM Single repeat found...");
						// PROBLEM
					}

					// calculate the score of wrap-around repeats
					// float repeatsScore = Trace.score(profile, sequence,
					// subst, 0, 0, repeats); // profile should never be scored
					// with gaps...
					float repeatsScore = repeats.size()
							* Trace.score(profile, sequence, subst, 0, 0,
									repeats); // profile should never be scored
												// with gaps...

					if (repeatsScore > bestRepeatsScore) {
						bestRepeatsScore = repeatsScore;
						bestRepeats = repeats;
						bestProfile = profile;
						bestEvdParams = evdParams;
						bestStartWrapping = startWrapping;
					}
				}

			} else {
				// no repeats detected: don't admit it
				// reportRepeats(new Sequence(""), 0, 0, new Vector(),
				// proteinName);
			}
		}

		if (bestRepeats.size() > 1) {
			Sequence sequenceOrig = new Sequence(seqOrig);
			// if (html)
			// reportRepeatsHtml(sequenceOrig, bestProfile.length(),
			// bestStartWrapping, bestRepeats, proteinName, typeNum);
			// else
			 reportMSA(sequenceOrig, bestProfile.length(),
					 bestStartWrapping, bestRepeats, proteinName, typeNum);
			 setRepeats_output(bestRepeats);

			if (DEBUG_VISUAL_SUMMARY) {
				System.err
						.println("Warning profile is not corrected with wrapa algorithm"); // TODO
				// System.err.println(bestProfile.toStringRepresentation());
			}

		} else {
			if (bestProfile != null && bestRepeats.size() == 0) {
				// System.err.println("Despite a profile was extracted, "
				// + bestRepeats.size()
				// + " repeats found using the profile " + proteinName);
			}
		}

		pb.stop();
		return bestRepeats;
	}

	private static void printTraces(String seq, Collection traces) {
		int traceno = 1;
		for (Iterator iter = traces.iterator(); iter.hasNext();) {
			Trace t = (Trace) iter.next();
			// System.err.println("Trace " + traceno);
			t.show(seq, seq, System.err);
			// System.err.println(t.toString(seq, seq));
			traceno++;
		}
	}

	/**
	 * @param startWrapping
	 * @param nonWrapa
	 * @return Repeats longer than half of the profile
	 */
	private static Collection filterTooShortRepeats(Collection allRepeats,
			int profileLen) {
		Collection longRepeats = new Vector();
		for (Iterator iter = allRepeats.iterator(); iter.hasNext();) {
			Trace t = (Trace) iter.next();
			if (t.length() > profileLen / 2)
				longRepeats.add(t);
		}
		return longRepeats;
	}

	/**
	 * @param startWrapping
	 * @param nonWrapa
	 * @return Repeats longer than half of the profile
	 */
	private static Collection filterTooWeakRepeats(Collection allRepeats,
			ComparableSequence cmpSeqX, ComparableSequence cmpSeqY,
			SubstitutionTable subst, float gapo, float gapx, double pvalue,
			double evdK, double evdLambda, String proteinName) {
		Collection longRepeats = new Vector();
		for (Iterator iter = allRepeats.iterator(); iter.hasNext();) {
			Trace t = (Trace) iter.next();

			// calculate teh
			float value = t.score(cmpSeqX, cmpSeqY, subst, gapo, gapx);
			float valueng = t.score(cmpSeqX, cmpSeqY, subst, 0, 0);

			// no gap is estimated statistically!
			double bychance = EstimateEVDParams.calculatePValue(1, valueng,
					cmpSeqX.length(), cmpSeqY.length(), false, evdK, evdLambda);
			// System.err.print("# For protein " + proteinName + " repeat "
			// + t.end[0] + ".." + t.end[t.end.length - 1] + " scored "
			// + value + " scored (no gaps) " + valueng);
			if (bychance <= pvalue) {
				// System.err.print(" included");
				longRepeats.add(t);
			} else {
				// System.err.print(" excluded");
			}
			// System.err.println();
		}
		return longRepeats;
	}

	private static void saveImageDecorated(float[][] m, String proteinName,
			String type, String outputDir) {

		float deco1[][] = Matrix.copy(m);
		Matrix.mirrorTopRight(deco1);
		float deco2[][] = Matrix.copy(deco1);
		Matrix.decorate(deco1, true);
		Matrix.decorate(deco2, true);
		// Matrix.decorateWithGrid(deco1);
		Matrix.saveImageRGB(deco2, deco2, deco1, VISUAL_PICT_INVERT,
				createImageFileName(outputDir, proteinName, type));
	}

	private static int calcIntRepeatsFromWrapa(Collection repeats) {
		return ((Trace) repeats.iterator().next()).start[0];
		/*
		 * TODO // determine the start of the repeat // possible starts of the
		 * repeat Collection possibleStarts = new Vector(); for (Iterator iter =
		 * repeats.iterator(); iter.hasNext();) { Trace t = (Trace) iter.next();
		 * 
		 * int start = t.start[0]; Integer startInteger = new Integer(start); if
		 * (!possibleStarts.contains(startInteger))
		 * possibleStarts.add(startInteger); }
		 * 
		 * for (Iterator iter = possibleStarts.iterator(); iter.hasNext();) {
		 * Integer startInteger = (Integer) iter.next(); int start =
		 * startInteger.intValue();
		 * 
		 * // calculate the integer number of repeats, // counting more than
		 * SelfSimilarity.PARAM_MIN_REPEAT of repeat // as an additional repeat
		 * calcIntRepeatsFromWrapa(start, repeats); }
		 */
	}


	private void reportMSA(ComparableSequence filteredSeq,
									  int profileLength, int _startWrapping, Collection repeats,
									  String proteinName, int typeNum) {

		// sanity check - repeats are to be extracted from the sequence
		Assert.assertTrue(filteredSeq.isSequence());

		// 0-based start of wrapping
		int startWrapping0 = _startWrapping - 1;

		// retrieve the sequence
		String seq = ((Sequence) filteredSeq).toString();

		// write down all the repeat instances
		if (repeats.size() == 0)
			return;
//		System.out.println("# type of the repeat");
//		System.out.println("REPEAT_TYPE " + typeNum);
//		System.out.println("# profile_length");
//		System.out.println("REPEAT_LENGTH " + profileLength);
//		System.out.println("# The list of repeats in the format:");
//		System.out.println("# START LENGTH [PVALUE [SCORE]]");
		List<String> msa = new ArrayList<String>();
		List<String> units = new ArrayList<String>();
		int repeatNo = 1;
		for (Iterator iter = repeats.iterator(); iter.hasNext(); repeatNo++) {
			Trace t = (Trace) iter.next();
			int traceMatchNum = t.end.length;
//			System.out.print(t.end[0] + " "
//					+ (t.end[traceMatchNum - 1] - t.end[0] + 1));
//			System.out.print("\t# Repeat " + repeatNo);
//			System.out.println();
			int start = t.end[0] - 1;
			int end = t.end[traceMatchNum - 1] -1 ;
			units.add(start+"-"+end);

		}

		// System.out.println("//");
//		System.out.println("# The multiple alignment of repeats");
//		System.out.println("# lo-case letters: not a part of alignment");
		// calculate profile gaps, 0-based
		int[] profGaps = Trace.getYGaps(profileLength, repeats);
//		Assert.assertTrue("no gap before the first position allowed",
//				profGaps[startWrapping0] == 0);
//		Assert.assertTrue(
//				"For wrapa traces, there should never be gap after the last pos",
//				profGaps[profGaps.length - 1] == 0);
		// print the "profile"
//		System.out
//				.println("# Profile pattern, \"X\": profile column, \"-\": a gap");
//		System.out.print("# ");
//		for (int i = startWrapping0;;) {
//			for (int j = 0; j < profGaps[i]; j++) {
//				System.out.print("-");
//			}
//			System.out.print("X");
//			i = Matrix.wrapNextPos(i, profileLength);
//			if (i == startWrapping0)
//				break;
//		}
//		System.out.println();

		// write down all the repeat instances

		// prepare gaps
		int gapMax = Matrix.getMax(profGaps) + profileLength;
		StringBuffer gapsAlign = new StringBuffer(gapMax);
		for (int j = 0; j < gapMax; j++) {
			gapsAlign.append("-");
		}
		// String gapsAlign = gapsAlignsb.toString();

		// iterate over the repeats
		repeatNo = 1;
		for (Iterator iter = repeats.iterator(); iter.hasNext(); repeatNo++) {
			Trace t = (Trace) iter.next();

			// fasta header
//			System.out.println(">Repeat " + repeatNo);

			String repeatAlign = createAlignment(t, seq, startWrapping0,
					profileLength, profGaps, gapsAlign);
			msa.add(repeatAlign);

//			System.out.println(repeatAlign);
		}

		RepeatSeq repeatSeq = new RepeatSeq();
		for (int i = 0; i < units.size() ; i++) {
			String unit = units.get(i);
			int start = Integer.parseInt(unit.split("-")[0]);
			int end = Integer.parseInt(unit.split("-")[1]);
			RepeatUnitSeq ru = new RepeatUnitSeq(msa.get(i),start,end);
			repeatSeq.getUnits().add(ru);
		}

		repeatSqs.add(repeatSeq);

		// System.out.println("//");
	}


	/*
	 * @param start
	 * 
	 * @param repeats
	 *//*
		 * private static int calcIntRepeatsFromWrapa(int start, Collection
		 * repeats) { for (Iterator iter = repeats.iterator(); iter.hasNext();)
		 * { Trace t = (Trace) iter.next();
		 * 
		 * 
		 * } }
		 */

	private static void reportRepeats(ComparableSequence filteredSeq,
			int profileLength, int _startWrapping, Collection repeats,
			String proteinName, int typeNum) {

		// sanity check - repeats are to be extracted from the sequence
		Assert.assertTrue(filteredSeq.isSequence());

		// 0-based start of wrapping
		int startWrapping0 = _startWrapping - 1;

		// retrieve the sequence
		String seq = ((Sequence) filteredSeq).toString();

		// write down all the repeat instances
		if (repeats.size() == 0)
			return;
		System.out.println("# type of the repeat");
		System.out.println("REPEAT_TYPE " + typeNum);
		System.out.println("# profile_length");
		System.out.println("REPEAT_LENGTH " + profileLength);
		System.out.println("# The list of repeats in the format:");
		System.out.println("# START LENGTH [PVALUE [SCORE]]");
		int repeatNo = 1;
		for (Iterator iter = repeats.iterator(); iter.hasNext(); repeatNo++) {
			Trace t = (Trace) iter.next();
			int traceMatchNum = t.end.length;
			System.out.print(t.end[0] + " "
					+ (t.end[traceMatchNum - 1] - t.end[0] + 1));
			System.out.print("\t# Repeat " + repeatNo);
			System.out.println();
		}

		// System.out.println("//");
		System.out.println("# The multiple alignment of repeats");
		System.out.println("# lo-case letters: not a part of alignment");
		// calculate profile gaps, 0-based
		int[] profGaps = Trace.getYGaps(profileLength, repeats);
		Assert.assertTrue("no gap before the first position allowed",
				profGaps[startWrapping0] == 0);
		Assert.assertTrue(
				"For wrapa traces, there should never be gap after the last pos",
				profGaps[profGaps.length - 1] == 0);
		// print the "profile"
		System.out
				.println("# Profile pattern, \"X\": profile column, \"-\": a gap");
		System.out.print("# ");
		for (int i = startWrapping0;;) {
			for (int j = 0; j < profGaps[i]; j++) {
				System.out.print("-");
			}
			System.out.print("X");
			i = Matrix.wrapNextPos(i, profileLength);
			if (i == startWrapping0)
				break;
		}
		System.out.println();

		// write down all the repeat instances

		// prepare gaps
		int gapMax = Matrix.getMax(profGaps) + profileLength;
		StringBuffer gapsAlign = new StringBuffer(gapMax);
		for (int j = 0; j < gapMax; j++) {
			gapsAlign.append("-");
		}
		// String gapsAlign = gapsAlignsb.toString();

		// iterate over the repeats
		repeatNo = 1;
		for (Iterator iter = repeats.iterator(); iter.hasNext(); repeatNo++) {
			Trace t = (Trace) iter.next();

			// fasta header
			System.out.println(">Repeat " + repeatNo);

			String repeatAlign = createAlignment(t, seq, startWrapping0,
					profileLength, profGaps, gapsAlign);

			System.out.println(repeatAlign);
		}
		// System.out.println("//");
	}

	private static void reportRepeatsHtml(ComparableSequence filteredSeq,
			int profileLength, int _startWrapping, Collection repeats,
			String proteinName, int typeNum) {

		// sanity check - repeats are to be extracted from the sequence
		Assert.assertTrue(filteredSeq.isSequence());

		// 0-based start of wrapping
		int startWrapping0 = _startWrapping - 1;

		// retrieve the sequence
		String seq = ((Sequence) filteredSeq).toString();

		// write down all the repeat instances
		if (repeats.size() == 0)
			return;

		System.out.println("<p><table border>");
		System.out.println("<tr class=\"darkgrey\">");
		System.out.println("<td colspan=\"4\" pairAlign=\"center\">Repeat type "
				+ typeNum + "</td>");
		System.out.println("</tr>");
		System.out.println("<tr class=\"lightgrey\">");
		System.out.println("<td>id</td>");
		System.out.println("<td>sequence</td>");
		System.out.println("<td>start</td>");
		System.out.println("<td>size</td>");
		System.out.println("</tr>");
		String aligns[] = repeatAlign(profileLength, repeats, startWrapping0,
				seq);
		int rstart[] = new int[repeats.size()];
		int rsize[] = new int[repeats.size()];

		int repeatNo = 1;
		for (Iterator iter = repeats.iterator(); iter.hasNext(); repeatNo++) {
			Trace t = (Trace) iter.next();
			int traceMatchNum = t.end.length;
			rstart[repeatNo - 1] = t.end[0];
			rsize[repeatNo - 1] = t.end[traceMatchNum - 1] - t.end[0] + 1;
		}

		repeatNo = 1;
		for (Iterator iter = repeats.iterator(); iter.hasNext(); repeatNo++) {
			Trace t = (Trace) iter.next();
			System.out.println("<tr>");
			System.out.println("<td>" + repeatNo + "</td>");
			System.out.println("<td><tt>" + aligns[repeatNo] + "</tt></td>");
			System.out.println("<td>" + rstart[repeatNo - 1] + "</td>");
			System.out.println("<td>" + rsize[repeatNo - 1] + "</td>");
			System.out.println("</tr>");
		}

		System.out.println("<table>");

		System.out.println("In the sequence:<br>");
		String fname = "repeatimg" + typeNum;
		writeRepeatImg(rstart, rsize, filteredSeq.length(), fname);

		System.out.println("<img src=\"" + fname + ".png\"><br>");
	}

	private static void writeRepeatImg(int[] rstart, int[] rsize,
			int protein_len, String fname) {
		final int width = 640;
		final int height = 20;
		final int fs = 2, fe = height - 3;

		float img_r[][] = new float[height][width], img_g[][] = new float[height][width], img_b[][] = new float[height][width];

		final float bg_r = 0.5f, bg_g = 0.3f, bg_b = .1f;
		final float re_r = .7f, re_g = .7f, re_b = 0;
		final float fi_r = .5f, fi_g = 0.5f, fi_b = 0.1f;

		Matrix.fill(img_g, bg_g);
		Matrix.fill(img_r, bg_r);
		Matrix.fill(img_b, bg_b);

		for (int i = 0; i < rstart.length; i++) {
			int s = (rstart[i] * width) / protein_len;
			int e = ((rstart[i] + rsize[i]) * width) / protein_len - 1;

			Matrix.drawVLine(img_r, s, fs, s, fe, re_r);
			Matrix.drawVLine(img_g, s, fs, s, fe, re_g);
			Matrix.drawVLine(img_b, s, fs, s, fe, re_b);
			Matrix.drawVLine(img_r, e, fs, e, fe, re_r);
			Matrix.drawVLine(img_g, e, fs, e, fe, re_g);
			Matrix.drawVLine(img_b, e, fs, e, fe, re_b);
			Matrix.drawHLine(img_r, s, fs, e, fs, re_r);
			Matrix.drawHLine(img_g, s, fs, e, fs, re_g);
			Matrix.drawHLine(img_g, s, fs, e, fs, re_g);
			Matrix.drawHLine(img_r, s, fe, e, fe, re_r);
			Matrix.drawHLine(img_g, s, fe, e, fe, re_g);
			Matrix.drawHLine(img_g, s, fe, e, fe, re_g);

			for (int x = s + 1; x < e; x++) {
				Matrix.drawVLine(img_r, x, fs + 1, x, fe - 1, fi_r);
				Matrix.drawVLine(img_g, x, fs + 1, x, fe - 1, fi_g);
				Matrix.drawVLine(img_b, x, fs + 1, x, fe - 1, fi_b);
			}

		}

		Matrix.saveImageRGBMinMax(img_r, img_g, img_b, 0f, 1f, false, fname);
	}

	/**
	 * @param profileLength
	 * @param repeats
	 * @param startWrapping0
	 * @param seq
	 */
	private static String[] repeatAlign(int profileLength, Collection repeats,
			int startWrapping0, String seq) {
		int repeatNo;
		String aligns[] = new String[repeats.size() + 1];

		// calculate profile gaps, 0-based
		int[] profGaps = Trace.getYGaps(profileLength, repeats);
		Assert.assertTrue("no gap before the first position allowed",
				profGaps[startWrapping0] == 0);
		Assert.assertTrue(
				"For wrapa traces, there should never be gap after the last pos",
				profGaps[profGaps.length - 1] == 0);
		// print the "profile"

		StringBuffer sb = new StringBuffer();

		for (int i = startWrapping0;;) {
			for (int j = 0; j < profGaps[i]; j++) {
				sb.append("-");
			}
			sb.append("X");
			i = Matrix.wrapNextPos(i, profileLength);
			if (i == startWrapping0)
				break;
		}
		aligns[0] = sb.toString();

		// write down all the repeat instances

		// prepare gaps
		int gapMax = Matrix.getMax(profGaps) + profileLength;
		StringBuffer gapsAlign = new StringBuffer(gapMax);
		for (int j = 0; j < gapMax; j++) {
			gapsAlign.append("-");
		}
		// String gapsAlign = gapsAlignsb.toString();

		// iterate over the repeats
		repeatNo = 1;
		for (Iterator iter = repeats.iterator(); iter.hasNext(); repeatNo++) {
			Trace t = (Trace) iter.next();

			String repeatAlign = createAlignment(t, seq, startWrapping0,
					profileLength, profGaps, gapsAlign);

			aligns[repeatNo] = repeatAlign;
		}
		return aligns;
	}

	public static String createAlignment(Trace t, String seq,
			int startWrapping0, int profileLength, int[] profGaps,
			StringBuffer gapsAlign) {
		Trace t0 = (Trace) t.clone();
		t0.move(-1, -1); // make the life easier by ommiting all exp-1
		return createAlignment0(t0, seq, startWrapping0, profileLength,
				profGaps, gapsAlign);
	}

	public static String createAlignment0(Trace t0, String seq,
			int startWrapping0, int profileLength, int[] profGaps,
			StringBuffer gapsAlign) {

		int tLen = t0.end.length;
		// write down the repeat instance
		StringBuffer repeatAlign = new StringBuffer(profileLength);
		// iterate through the profile
		for (int i = startWrapping0, tPos = 0;;) {

			// number of gaps according to the profile
			int posGapsProf = profGaps[i];

			// calulate the gap size in seq
			int posGapsSeq = 0;
			if (tPos > 0 && tPos < tLen) {
				posGapsSeq = (t0.end[tPos] - t0.end[tPos - 1]) - 1;
				// insert skipped residues from seq
				repeatAlign.append(seq.substring(t0.end[tPos - 1] + 1,
						t0.end[tPos - 1] + 1 + posGapsSeq).toLowerCase());
			}
			// fill in the rest with gaps
			repeatAlign
					.append(gapsAlign.substring(0, posGapsProf - posGapsSeq));
			// sanity check
			if (Assert.DEBUG) {
				if (posGapsSeq > 0)
					Assert.assertTrue(i == t0.start[tPos]);
			} // make sure that when there was a gap in the seq, the letter will
				// be printed

			// if this is a position to print a letter
			if (tPos < tLen && i == t0.start[tPos]) {
				// print the letter
				repeatAlign.append(Character.toUpperCase(seq
						.charAt(t0.end[tPos])));
				tPos++;
			} else {
				// there is a gap
				repeatAlign.append("-");
			}

			// next iteration
			int nexti = Matrix.wrapNextPos(i, profileLength);
			if (nexti == startWrapping0)
				break;
			i = nexti;
		}

		if (Assert.DEBUG) {
			int gapSum = 0;
			for (int j = 0; j < profGaps.length; j++) {
				gapSum += profGaps[j];
			}
			Assert.assertTrue("Improper pairAlign length",
					repeatAlign.length() == profileLength + gapSum);
		}

		if (Assert.DEBUG) {
			// double check the alignment contents
			int tPos = 0;
			int aPos = 0;
			while (aPos < repeatAlign.length()) {
				char c = repeatAlign.charAt(aPos);
				if (Character.isUpperCase(c)) {
					Assert.assertTrue(Character.toUpperCase(seq
							.charAt(t0.end[tPos])) == c);
					tPos++;
				}
				aPos++;
			}
			Assert.assertTrue(tPos == t0.length());
		}

		return repeatAlign.toString();
	}

	/**
	 * @see Extract method
	 */
	public static class RepeatStatisticsParam {

		// double maxSWCoverageRatio;
		double maxAllTraceLen;

		public RepeatStatisticsParam(int maxAllTraceLen/*
														 * , double
														 * maxSWCoverageRatio
														 */) {
			this.maxAllTraceLen = maxAllTraceLen;
			// this.maxSWCoverageRatio = maxSWCoverageRatio;
		}

	}

	/**
	 * Calculate maximal tandem repeat size based on the trace with maximal
	 * length.
	 * 
	 * @param
	 * @return The area which the strongest trace covers
	 */
	/*
	 * private static int calculateMaxTandemRepeatSize(Collection bestTraces) {
	 * int maxTandemRepeat = 0; for (Iterator iter = bestTraces.iterator();
	 * iter.hasNext();) { Trace element = (Trace) iter.next(); maxTandemRepeat =
	 * Math.max(maxTandemRepeat, element.covered()); //!!! to fix - it's not
	 * really the strongest } return maxTandemRepeat; }
	 */

	private static Collection/* int */calculateTandemRepeatSize(
			float[][] transM, int maxTandemRepeatSize, int start, int end,
			String outputDir, String proteinName, int typeNum) {
		// can happen: Assert.assertTrue(end-start+1 >= maxTandemRepeatSize);
		// -- the whole trace is projected onto single point
		float[] _projection = Matrix.projectTriangleSE(transM, start, end
				- start + 1);

		Vector rSizes = new Vector();

		// smooth the projection accordingly to the param
		float[] projectionAvg = Matrix.createMovingAverage(_projection,
				PARAM_SMOOTH_PROJECTION);

		if (VISUAL_PROJECT_SE && DEBUG_VISUAL_SUMMARY) {
			float[][] histM = Matrix.paintHistogram(Matrix.createMovingAverage(
					Matrix.projectTriangleSE(transM), PARAM_SMOOTH_PROJECTION),
					null, 100);
			Matrix.saveImageGray(
					histM,
					VISUAL_PICT_INVERT,
					createImageFileName(outputDir, proteinName, "SW proj "
							+ typeNum));
		}

		// take the strongest trace peak in SE projection
		int maxTraceOffset = Matrix.getMaxPos(projectionAvg); // distance to the
																// diagonal
		rSizes.add(0, new Integer(maxTraceOffset));
		int smoothMaxTraceOffset = maxTraceOffset;
		float maxTrace = projectionAvg[maxTraceOffset];

		// move to the left with 1/4-3/4 method
		float tmpMax = maxTrace / 2;
		int lowerBound = maxTraceOffset / 4;
		int upperBound = ((maxTraceOffset + 1) * 3) / 4; // PROBLEM +1 goes
															// outside!
															// ((maxTraceOffset)*3)/4
															// + 1
		upperBound = Math.max(upperBound, maxTraceOffset - 1);
		for (int i = lowerBound; i <= upperBound; i++) {
			if (projectionAvg[i] >= tmpMax && maxTraceOffset != i) {
				if (DEBUG_VISUAL_SIMPLE) {
					// System.err.println("Max offset moved from "
					// + maxTraceOffset + " to " + i);
				}
				rSizes.add(0, new Integer(i));
				maxTraceOffset = i;
				tmpMax = projectionAvg[i];
				break;
			}
		}
		// Assert.assertTrue(Math.min(maxTraceOffset+1, maxTandemRepeatSize) ==
		// ((Integer) rSizes.elementAt(0)).intValue());

		// move right from the current highest
		maxTraceOffset = ((Integer) rSizes.elementAt(0)).intValue();
		tmpMax = projectionAvg[maxTraceOffset] / 2;

		lowerBound = (maxTraceOffset * 3) / 2 + 1;
		if (lowerBound <= projectionAvg.length - 1) {
			maxTraceOffset = Matrix.getMaxPos(projectionAvg, lowerBound,
					projectionAvg.length - 1);

			if (projectionAvg[maxTraceOffset] >= tmpMax
					&& maxTraceOffset <= maxTandemRepeatSize) {
				rSizes.add(new Integer(maxTraceOffset));
			}
		}

		/*
		 * // retrieve the exact maxTraceOffset, which is in neighborhood of max
		 * int oldMaxTraceOffset = maxTraceOffset; lowerBound =
		 * maxTraceOffset-PARAM_SMOOTH_PROJECTION/2; lowerBound =
		 * Math.max(lowerBound, 0); upperBound =
		 * maxTraceOffset-PARAM_SMOOTH_PROJECTION/2+PARAM_SMOOTH_PROJECTION-1;
		 * //!!!maxTraceOffset = Matrix.getMaxPos(_projection, lowerBound,
		 * upperBound); if (DEBU_G_REPEAT_LEN_RETRIVAL && maxTraceOffset !=
		 * oldMaxTraceOffset)
		 * //System.err.println("Max offset moved from "+oldMaxTraceOffset
		 * +" to "+maxTraceOffset);
		 */
		// tandem repeat size is the max found +1
		for (int i = 0; i < rSizes.size(); i++) {
			int v = ((Integer) rSizes.elementAt(i)).intValue();
			v = Math.min(v + 1, maxTandemRepeatSize);
			rSizes.set(i, new Integer(v));
		}
		// filter the values, leaving only unique ones
		for (int i = 0; i < rSizes.size(); i++) {
			Integer v = (Integer) rSizes.elementAt(i);
			int j;
			while ((j = rSizes.indexOf(v, i + 1)) != -1) {
				rSizes.remove(j);
			}
		}
		// tandem repeat size is the max found +1
		// return Math.min(maxTraceOffset+1, maxTandemRepeatSize);
		return rSizes;
	}

	/**
	 * @param transM
	 * @param startx
	 * @param endx
	 * @param starty
	 * @param endy
	 * @param maxTraceOffset
	 * @return Max scoring region
	 */
	private static int maxScoringStrap(float[][] m, int start, int end,
			int tandemRepeatSize, int distanceFromDiagonal) {
		if (tandemRepeatSize == 0)
			return start;
		// find maximal noise level
		float maxNoise = 0;
		/*
		 * for (int y = start; y < end; y++) { maxNoise = Math.max(maxNoise,
		 * sumNoise(m, y, tandemRepeatSize)); }
		 */
		float sum = 0;
		// calculate initial max strap
		for (int y = start; y < start + tandemRepeatSize; y++) {
			sum += sumLine(m, y, start, end, tandemRepeatSize,
					distanceFromDiagonal, maxNoise);
		}
		float maxSum = sum;
		int maxSumPos = start;
		// advance the strap down
		for (int y = start; y <= end - tandemRepeatSize; y++) {
			sum -= sumLine(m, y, start, end, tandemRepeatSize,
					distanceFromDiagonal, maxNoise);
			sum += sumLine(m, y + tandemRepeatSize, start, end,
					tandemRepeatSize, distanceFromDiagonal, maxNoise);
			if (sum > maxSum) {
				maxSum = sum;
				maxSumPos = y + 1;
			}
		}
		return maxSumPos;
	}

	/**
	 * Sums the horizontal line y, starting from x=start till end
	 * 
	 * @param m
	 * @param y
	 * @param start
	 * @param end
	 * @return Calculated sum
	 */
	private static float sumLine(float[][] m, int y, int start, int end,
			int tandemRepeatSize, int distanceFromDiagonal, float maxNoise) {
		// the code below is mainly a copy of Profile.Profile
		// don't forget to update both if you change something!

		Assert.assertTrue("Matrix not mirrored", Matrix.isMirrored(m));
		float sum = 0;

		int maxSize = Profile.getMaxPosWithFilteringArraySize(start, end,
				tandemRepeatSize);
		int[] pos = new int[maxSize];
		float[] val = new float[maxSize];

		int nPos = Profile.getMaxPosWithFiltering(m, y, start, end,
				tandemRepeatSize, distanceFromDiagonal, pos, val, null);

		for (int i = 0; i < nPos; i++) {
			sum += val[i];
		}

		return sum;
	}

	/**
	 * Sums the horizontal line y, starting from x=start till end
	 * 
	 * @param m
	 * @param y
	 * @param start
	 * @param end
	 * @return Calculated sum
	 */
	private static float _sumLine(float[][] m, int y, int start, int end,
			int tandemRepeatSize, float maxNoise) {
		Assert.assertTrue("Matrix not mirrored", Matrix.isMirrored(m));
		float sum = 0;

		for (int x = start; x <= end; x++) {
			sum += m[y][x];
		}

		// scale the value
		float sumWideTrace;
		float sumNoise;
		int traceMargin;

		// f(1) = 0;
		// f(2) = 0;
		// f(3) = 0;
		// f(4) = 0;
		// f(5) = 1;
		// f(6) = 1;
		// f(7) = 1;
		// f(8) = 2;
		// f(9) = 2;
		// traceMargin = calcTraceMargin(tandemRepeatSize);

		// sumWideTrace = sumWideTrace(m, y, tandemRepeatSize);
		float noiseLevel;
		sumNoise = sumNoise(m, y, tandemRepeatSize);

		if (maxNoise != 0f)
			noiseLevel = sumNoise / maxNoise;
		else
			noiseLevel = 0.0f;
		noiseLevel = Math.max(noiseLevel, 0.0f);
		noiseLevel = Math.min(noiseLevel, 0.5f);
		// //System.err.println(noiseLevel);

		return sum * (1.0f - noiseLevel);
	}

	private static int calcTraceMargin(int tandemRepeatSize) {
		int traceMargin;
		traceMargin = (tandemRepeatSize + 1) / PARAM_WIDE_TRACE - 1;
		traceMargin = Math.max(traceMargin, 0);
		Assert.assertTrue(traceMargin < tandemRepeatSize);
		return traceMargin;
	}

	private static float sumWideTrace(float[][] m, int y, int tandemRepeatSize) {

		float sumWideTrace;
		sumWideTrace = 0f;
		int traceMargin = calcTraceMargin(tandemRepeatSize);
		for (int x = tandemRepeatSize - traceMargin; x <= tandemRepeatSize
				+ traceMargin; x++) {
			Assert.assertTrue(x > 0);
			Assert.assertTrue(x < tandemRepeatSize * 2);
			sumWideTrace += m[y][y + x];
		}
		return sumWideTrace;
	}

	/**
	 * Sums the noise level from diagonal to the trandemRepeatSize-traceMargin,
	 * in both directions.
	 * 
	 * @param m
	 * @param y
	 * @param tandemRepeatSize
	 * @return
	 */
	private static float sumNoise(float[][] m, int y, int tandemRepeatSize) {

		Assert.fail("Shouldn't be used!");
		float sumNoise;
		sumNoise = 0f;
		int traceMargin = calcTraceMargin(tandemRepeatSize);
		for (int x = -(tandemRepeatSize - traceMargin) + 1; x < tandemRepeatSize
				- traceMargin; x++) {
			// Assert.assertTrue(x > 0);
			Assert.assertTrue(x < tandemRepeatSize);
			Assert.assertTrue(x > -tandemRepeatSize);
			if (y + x >= 0 && y + x < m[0].length)
				sumNoise += m[y][y + x];
		}

		return sumNoise;
	}

	private static void createHistograms(float[][] m, String proteinName,
			String suffix, String outputDir) {
		float[][] histM;
		if (VISUAL_PROJECT_E) {
			histM = Matrix.paintHistogram(Matrix.projectE(m), null, 100);
			Matrix.saveImageGray(
					histM,
					VISUAL_PICT_INVERT,
					createImageFileName(outputDir, proteinName, "histE "
							+ suffix));
		}
		if (VISUAL_PROJECT_SW) {
			histM = Matrix.paintHistogram(projectSW(m), null, 100);
			Matrix.saveImageGray(
					histM,
					VISUAL_PICT_INVERT,
					createImageFileName(outputDir, proteinName, "histSW "
							+ suffix));
		}
		if (VISUAL_PROJECT_SE) {
			histM = Matrix.paintHistogram(Matrix.projectTriangleSE(m), null,
					100);
			Matrix.saveImageGray(
					histM,
					VISUAL_PICT_INVERT,
					createImageFileName(outputDir, proteinName, "histSE "
							+ suffix));

			histM = Matrix.paintHistogram(Matrix.createMovingAverage(
					Matrix.projectTriangleSE(m), PARAM_SMOOTH_PROJECTION),
					null, 100);
			Matrix.saveImageGray(
					histM,
					VISUAL_PICT_INVERT,
					createImageFileName(outputDir, proteinName, "histSE avg "
							+ PARAM_SMOOTH_PROJECTION + " " + suffix));
		}
	}

	/**
	 * @param bestTraces
	 * @return
	 */
	private static Vector convertToBinary(Vector bestTraces) {
		Vector binaryTraces = new Vector(bestTraces.size());
		for (int i = 0; i < bestTraces.size(); i++) {
			Trace t = (Trace) ((Trace) bestTraces.elementAt(i)).clone();
			t.setValue(1.0f);
			binaryTraces.add(t);
		}
		return binaryTraces;
	}

	/**
	 * Recalculates the weights of trace taking aa similarity into account.
	 * 
	 * @param bestTraces
	 * @return
	 */
	private static Vector convertToSensitive(Vector bestTraces, AlignData data) {
		Vector sensitiveTraces = new Vector(bestTraces.size());
		// score = normalizedScore + local strength
		for (int i = 0; i < bestTraces.size(); i++) {
			Trace t = (Trace) ((Trace) bestTraces.elementAt(i)).clone();
			float score = t.score(data.getSeqX(), data.getSeqY(),
					data.getSubst(), data.getGapO(), data.getGapE());
			// make sensitive
			t.convertToSensitive(PARAM_SENSITIVE_TRACE_SMOOTH, data.getSeqX(),
					data.getSeqX(), data.getSubst(), data.getGapO(),
					data.getGapE());
			// normalize
			Matrix.correctScaleToSum(t.value, score);
			t.addValue((float) (score / Math.sqrt(t.length()))); // PROBLEM
			// store
			sensitiveTraces.add(t);
		}
		return sensitiveTraces;
	}

	/**
	 * Projects top right half of the matrix (including the diagonal) onto a
	 * vector in SW direction
	 * 
	 * @return Projected vector in SW dir
	 */
	private static float[] projectSW(float[][] m) {
		float[] dproject = new float[m.length + m[0].length - 1];
		for (int x = 0; x < m.length + m[0].length; x++) {
			for (int i = 0; i <= x && i < m.length; i++) {
				if (x - i < m[0].length - 1)
					dproject[x] += m[i][x - i] + m[i][x - i + 1];
			}
		}
		return dproject;
	}

	/**
	 * Prints histogram onto the output
	 */
	private static void printHistogram(float[] hist) {
		for (int i = 0; i < hist.length; i++) {
			System.out.println(hist[i]);
		}
	}

	public static String createImageFileName(String outputDir, String fname,
			String suffix) {
		return outputDir + fname + " " + suffix;
	}

	/**
	 * Method selectSubtrace selects the subtrace of trace of length dist
	 * 
	 * @param trace
	 * @param dist
	 * @return int[][]
	 */
	static private int[][] selectSubtrace(int[][] trace, int dist) {
		return null;
	}

	/**
	 * Calculates the distance of trace to diagonal. In future !!!the distance
	 * may be weighted with similarity matrix.
	 * 
	 * @param data
	 * @param trace
	 * @return int
	 */
	static private int distanceFromDiagonal(/* float [][]m, */int[][] trace) {
		// calculate the avg distance from diagonal
		int sumdist = 0;
		for (int i = 0; i < trace[0].length; i++) {
			sumdist += trace[0][i] - trace[1][i]; // x_i - y_i
		}
		return sumdist / trace[0].length;
	}

	/**
	 * Shows Waterman-Eggert traces (more significant than specified
	 * <code>pscore</code>).
	 */
	public static void showRepeats(String seq, SubstitutionTable subst,
			int gapo, int gapx, float pscore, String fname) {
		// final int gap = -4;
		// AlignData data = new AlignData(seq, seq, gap, gap, subst);
		AlignDataAffineGotoh data = new AlignDataAffineGotoh(seq, seq, gapo,
				gapx, subst);
		LocalInit.localInit(data);
		SelfLocalFill.selfLocalFill(data);
		float source_matrix[][] = Matrix.copy(data.getMatrix());
		float traces[][] = Matrix.createEmpty(source_matrix);

		// get sequential max regions
		float value = pscore;
		for (int i = 0; i < 30 && value >= pscore; i++) {
			// reconstruct the array
			value = LocalAlignAffineGotoh.localAlignAffine(data);
			if (value > pscore) {
				// leave the trace
				data.markMaxRegion(traces, (int[][]) null, value);
				// prepare source array for the next step
				data.markMaxRegion(source_matrix, (int[][]) null, -Matrix.INF);
				data.setMatrixFrom(source_matrix);
			}
		}
		Matrix.saveImageGray(traces, true, fname + " traces");

		Matrix.mirrorTopRight(traces);
		Matrix.saveImageGray(makeTransitive(traces), true, fname
				+ " traces trans");
	}

	public static void showRepeatsInOnePass(String seq,
			SubstitutionTable subst, String fname, String suffix, int gapo,
			int gapx, boolean makeTransitive, Align align) {
		AlignData data = new AlignDataAffineGotoh(seq, seq, gapo, gapx, subst);
		// AlignData data = new AlignData(seq, seq, gapo, gapx, subst);
		align.align(data);

		float matrix[][] = data.getMatrix();
		Matrix.scale(matrix, 0f, 1.0f);
		Matrix.saveImageGrayMinMax(matrix, 0f, 1f, true, fname + " " + gapo
				+ " " + gapx + " " + suffix);
		if (makeTransitive) {
			Matrix.mirrorTopRight(matrix);
			float[][] transitive = makeTransitive(matrix);
			// Matrix.subtract(transitive, matrix);
			// Matrix.saveImageGrayMinMax(transitive, 0f,
			// Matrix.getMax(transitive), true, fname+" trans "+gap);
			Matrix.scale(transitive, 0f, 1.0f);
			Matrix.saveImageRGBMinMax(transitive, matrix,
					Matrix.createEmpty(matrix), 0f, 1f, true, fname + " "
							+ gapo + " " + gapx + " trans" + suffix);
			Matrix.saveImageGrayMinMax(transitive, 0f, 1, true, fname + " "
					+ gapo + " " + gapx + " trans" + " GRAY" + suffix);
			// Matrix.saveImageRGBMinMax(transitive, Matrix.create(matrix),
			// matrix, 0f, 1f, true, fname+" trans "+gap);
		}
	}

	/**
	 * Graphically shows valeys from S-W algorithm, and calculates transitive
	 * matrix.
	 */
	public static void showRepeatsInOnePass(String seq,
			SubstitutionTable subst, String fname, int gapo, int gapx,
			boolean makeTransitive) {
		// showRepeatsInOnePass(seq, subst, fname, gapo, gapx, makeTransitive,
		// new CompleteAlignSelfAffine());
		showRepeatsInOnePass(seq, subst, fname, "min", gapo, gapx,
				makeTransitive, new AlignDDPMin(new CompleteAlignSelfAffine()));

	}

	/**
	 * Graphically shows valeys from S-W algorithm, and calculates transitive
	 * matrix.
	 */
	public static void showRepeatsInOnePass(String seq,
			SubstitutionTable subst, String fname, boolean makeTransitive) {
		showRepeatsInOnePass(seq, subst, fname, -10, -2, false);
	}

	private static float[][] makeTransitive(float[][] m) {
		// O(n^3)
		float[][] transitive = Matrix.createEmpty(m);
		Assert.assertTrue(m.length == m[0].length);
		int len = m.length;
		// ProgressBar pb = new ProgressBar();
		// pb.start("Calculating transitivity", len);
		for (int y = 0; y < len; y++) {
			for (int x = y; x < len; x++) {
				// transitive[y][x] = m[y][x];
				for (int z = 0; z < len; z++) {
					// transitive[y][x] += m[y][z] * m[z][x];
					// if (z != x && z != y)
					// transitive[y][x] = Math.max(transitive[y][x],
					// (m[y][z]+m[z][x])/2);
					// transitive[y][x] = Math.max(transitive[y][x],
					// (float)Math.sqrt((m[y][z]*m[y][z]+m[z][x]*m[z][x])/2));
					// calculate
					// transitive[y][x] = Math.max(transitive[y][x],
					// Math.min(m[y][z], m[z][x]));
					float old = transitive[y][x];
					float e1 = m[y][z];
					float e2 = m[z][x];
					if (e1 > 0 && e2 > 0) {
						float emin = e1 < e2 ? e1 : e2;
						transitive[y][x] = emin > old ? emin : old;
					}
				}
			}
			// pb.advance();
		}
		// pb.stop();
		return transitive;
	}

	public static float calcRepeats(AlignData data, Align align) {
		LocalInit.localInit(data);
		SelfLocalFill.selfLocalFill(data);
		return align.align(data);
	}

	public static void showEVDDiscrepancy(ProteinLib library, String name,
			SubstitutionTable subst, int gapo, int gapx, String outputDir) {
		ProgressBar pb = new FakeProgressBar();
		pb.start("Proteins analyzed", 1);

		// Waterman-Eggert
		pb.advance(name);
		Assert.assertTrue(library.get(name) != null);
		EstimateEVDParams.findEVDDiscrepancy(library.get(name), subst, gapo,
				gapx, 0.005f, true, name, outputDir);

		pb.stop();
	}

	public static void showAllEVDDiscrepancy(ProteinLib library,
			SubstitutionTable subst, int gapo, int gapx, String outputDir,
			int procnr, int proctotal) {
		List proteins = new Vector(library.getAllProteins());
		String fname = "data/seqs350-1000all_anal.txt";
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(
					new FileInputStream(fname)));

			String pname = null;
			ProgressBar pb = DEBUG_VISUAL_SUMMARY ? ((ProgressBar) new TextProgressBar())
					: ((ProgressBar) new FakeProgressBar());
			pb.start("Protein analyzed", proteins.size());
			while ((pname = br.readLine()) != null) {
				if (pname.length() == 0 || pname.startsWith("#"))
					continue;

				// pb.info("To analyze: "+pname);
				showEVDDiscrepancy(library, pname, subst, gapo, gapx, outputDir);

			}
			pb.stop();
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}

	}

	public void showTransitivity(ProteinLib library, ProteinLib librarySeg,
			String name, int maxLen, SubstitutionTable subst, int gapo,
			int gapx, String outputDir, boolean html) {
		ProgressBar pb = new FakeProgressBar();
		pb.start("Proteins analyzed", 1);

		// Waterman-Eggert
		pb.advance(name);
		Assert.assertTrue(library.get(name) != null);
		Collection repeatsFound = null;
		String seqOrig = library.get(name);
		String seq = librarySeg.get(name);

		Assert.assertTrue(seq != null && seqOrig != null);

		if (seq.length() <= maxLen) {

			if (html)
				System.out.println("<p><font size=\"+2\">Protein " + "<b>"
						+ name + "</b></font>");
			else {
				// System.out.println(">" + name);
			}

			if (!subst.isAllowedSeq(seq)) {
				if (html)
					System.out
							.println("<p><font color=\"red\">The protein contains characters which are not amino acids</font>");
				else
					System.out
							.println("# the protein contains characters which are not amino acids");
			} else if (seq.length() > 0) {
				boolean anyRepeatFound = false;

				for (int typeNum = 1; repeatsFound == null
						|| repeatsFound.size() > 0; typeNum++) {

					// TODO because X regions can be reported as repeats
					// this method can report the same residue involved in
					// different repeats
					repeatsFound = findRepeats(seq, seqOrig, subst, gapo, gapx,
							0.01f, name, outputDir, typeNum, html);

					if (repeatsFound != null)
						anyRepeatFound = true;

					// mask the repeats from the previous iteration
					seq = Trace.maskSequence(seq, repeatsFound);
				}

				if (!anyRepeatFound)
					if (html)
						System.out
								.println("<p><font size=\"+2\">No repeats found</font>");
					else
						System.out.println("# no repeats found");

			}
			if (html) {
				System.out
						.println("<p><font style=\"font-family: Arial,Helvetica,sans-serif; font-size: 11px\">If you are experiencing problems with the server, please contact <b>radek [at] cs vu nl</b>.</font>");
			} else {
				// System.out.println("# end of the protein " + name);
				// System.out.println("//");
			}
		} else {
			// System.err.println("Protein " + name + " longer (" + seq.length()
			// + ")than max " + maxLen);
		}

		pb.stop();
	}

	public void showAllTransitivity(ProteinLib library, ProteinLib librarySeg,
			int maxLen, SubstitutionTable subst, int gapo, int gapx,
			String outputDir, int procnr, int proctotal, boolean html) {
		List proteins = new Vector(library.getAllProteins());

		// cut proteins to the desired size, depending on the number of
		// processes
		int startprotein = (procnr * proteins.size()) / proctotal;
		int endprotein = ((procnr + 1) * proteins.size()) / proctotal; // protein
																		// which
																		// is
																		// not
																		// supposed
																		// to be
																		// calculated

		proteins = proteins.subList(startprotein, endprotein);

		// sort the list
		Collections.sort(proteins, new LengthComparator());

		ProgressBar pb = DEBUG_VISUAL_SUMMARY ? ((ProgressBar) new TextProgressBar())
				: ((ProgressBar) new FakeProgressBar());
		pb.start("Protein analyzed", proteins.size());
		for (Iterator iter = proteins.iterator(); iter.hasNext();) {
			String seq = (String) iter.next();
			String name = (String) library.getFirstName(seq);
			pb.info("To analyze: " + name);
			showTransitivity(library, librarySeg, name, maxLen, subst, gapo,
					gapx, outputDir, html);
			pb.advance(name);
		}
		pb.stop();
	}

	static final String PARAM_LIB = "-fasta";
	static final String PARAM_MATRIX = "-matrix";
	static final String PARAM_OUTPUT_DIR = "-o";
	static final String PARAM_PROTEIN_NAME = "-protein";
	static final String PARAM_PROC_NR = "-procNr";
	static final String PARAM_PROC_TOTAL = "-procTotal";
	static final String PARAM_EVD_DISCREPANCY = "-evd";
	static final String PARAM_NO_SEG = "-noseg";
	static final String PARAM_GAPO = "-gapo";
	static final String PARAM_GAPX = "-gapx";
	static final String PARAM_MAX_LEN = "-max";
	static final String ARG_FORCE = "-force";
	static final String ARG_NO_FORCE = "-noforce";
	static final String ARG_HTML = "-html";

	private int getParamPos(String[] argv, String param) {
		for (int i = argv.length - 1; i >= 0; i--) {
			if (argv[i].equalsIgnoreCase(param))
				return i;
		}
		return -1;
	}

	public Collection getRepeatsOutput() {
		return repeats_output;
	}

	public void setRepeats_output(Collection repeats_output) {
		this.repeats_output = repeats_output;
	}

}
