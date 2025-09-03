package edu.put.ma;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import junitx.framework.FileAssert;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.Maps;

import edu.put.ma.archiver.Archiver;
import edu.put.ma.archiver.ArchiverFactory;
import edu.put.ma.model.input.AnalysisInputModel;
import edu.put.ma.model.input.AnalysisInputModelImpl;
import edu.put.ma.model.input.ComparisonInputModelImpl;
import edu.put.ma.model.input.DeformationProfileInputModelImpl;

public class AssessmentTest {

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private Assessment assessment;

    @Before
    public void setUp() {
        assessment = new Assessment();
    }

    @Test
    public void testPdbValidationOfSingleModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath("SINGLE-STRUCTURE-ANALYSIS/PDB-VALIDATION/models/model.pdb")
                .outputFilePath("validation.xml").command(CommandEnum.PDB_VALIDATION);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/PDB-VALIDATION/expected");
    }

    @Test
    public void testPdbValidationOfMultipleModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath("SINGLE-STRUCTURE-ANALYSIS/PDB-VALIDATION/models")
                .outputFilePath("validation.xml").command(CommandEnum.PDB_VALIDATION);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/PDB-VALIDATION/expected");
    }

    @Test
    public void testClashScoreOfSingleContinousModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath(
                        "SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/continuous-models/4_adamiak_1_rpr.pdb")
                .outputFilePath("4_adamiak_1_rpr.xml").considerAtomsSupportedByRNAPuzzlesOnly(true)
                .command(CommandEnum.Measure.CS);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/expected/continuous-models");
    }

    @Test
    public void testClashScoreOfMultipleContinousModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath(
                        "SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/continuous-models")
                .outputFilePath("4_adamiak_1_rpr.xml").considerAtomsSupportedByRNAPuzzlesOnly(true)
                .command(CommandEnum.Measure.CS);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/expected/continuous-models");
    }

    @Test
    public void testClashScoreOfSingleIncontinousModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath(
                        "SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/incontinuous-models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.xml").considerAtomsSupportedByRNAPuzzlesOnly(true)
                .command(CommandEnum.Measure.CS);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/expected/incontinuous-models");
    }

    @Test
    public void testClashScoreOfMultipleIncontinousModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath(
                        "SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/incontinuous-models")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.xml").considerAtomsSupportedByRNAPuzzlesOnly(true)
                .command(CommandEnum.Measure.CS);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SCORES/CLASH-SCORE/expected/incontinuous-models");
    }

    @Test
    public void testWholeSequenceOfSingleContinousModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath(
                        "SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/continuous-models/models/4_adamiak_1_rpr.pdb")
                .outputFilePath("4_adamiak_1_rpr.xml").command(CommandEnum.Sequence.SEQUENCE);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/expected/continuous-models");
    }

    @Test
    public void testWholeSequenceOfMultipleContinousModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath(
                        "SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/continuous-models/models")
                .outputFilePath("4_adamiak_1_rpr.xml").command(CommandEnum.Sequence.SEQUENCE);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/expected/continuous-models");
    }

    @Test
    public void testWholeSequenceOfSingleIncontinousModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath(
                        "SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.xml").command(CommandEnum.Sequence.SEQUENCE);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/expected/incontinuous-models");
    }

    @Test
    public void testWholeSequenceOfMultipleIncontinousModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath(
                        "SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/incontinuous-models/models")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.xml").command(CommandEnum.Sequence.SEQUENCE);
        analyse(analysisInputModelBuilder,
                "./SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/SEQUENCE/expected/incontinuous-models");
    }

    @Test
    public void testFragmentsOfSequenceOfSingleIncontinousModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath(
                        "SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/FRAGMENT/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.xml").command(CommandEnum.Sequence.FRAGMENT);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/FRAGMENT/expected");
    }

    @Test
    public void testFragmentsOfSequenceOfMultipleIncontinousModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath("SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/FRAGMENT/models")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.xml").command(CommandEnum.Sequence.FRAGMENT);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/SEQUENCE/FRAGMENT/expected");
    }

    @Test
    public void testOriginal3dOfSingleIncontinousModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath(
                        "SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/ORIGINAL-3D/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.zip").command(CommandEnum.Structure3d.ORIGINAL_3D);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/ORIGINAL-3D/expected");
    }

    @Test
    public void testOriginal3dOfMultipleIncontinousModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath("SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/ORIGINAL-3D/models")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.zip").command(CommandEnum.Structure3d.ORIGINAL_3D);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/ORIGINAL-3D/expected");
    }

    @Test
    public void testRenumerated3dOfSingleIncontinousModel() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .singlePdbModelFilePath(
                        "SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/RENUMERATED-3D/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.zip").command(CommandEnum.Structure3d.RENUMERATED_3D);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/RENUMERATED-3D/expected");
    }

    @Test
    public void testRenumerated3dOfMultipleIncontinousModels() throws Exception {
        final AnalysisInputModelImpl.Builder analysisInputModelBuilder = new AnalysisInputModelImpl.Builder();
        analysisInputModelBuilder
                .multiplePdbModelsDirectoryPath(
                        "SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/RENUMERATED-3D/models")
                .targetModelsAlignment("14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .outputFilePath("14_ChenPostExp_1_rpr.zip").command(CommandEnum.Structure3d.RENUMERATED_3D);
        analyse(analysisInputModelBuilder, "./SINGLE-STRUCTURE-ANALYSIS/STRUCTURE-3D/RENUMERATED-3D/expected");
    }

    @Test
    public void testRootMeanSquareDeviationOfSingleContinousModelWithoutAFewAtoms() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/lack-of-atoms/1_solution_0_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/lack-of-atoms/models/1_das_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("1_das_1_rpr.xml")
                .command(CommandEnum.Measure.RMSD);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/expected/lack-of-atoms");
    }

    @Test
    public void testRootMeanSquareDeviationOfMultipleContinousModelsWithoutAFewAtoms() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/lack-of-atoms/1_solution_0_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/lack-of-atoms/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("1_das_1_rpr.xml")
                .command(CommandEnum.Measure.RMSD);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/expected/lack-of-atoms");
    }

    @Test
    public void testRootMeanSquareDeviationOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .command(CommandEnum.Measure.RMSD);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/expected/continuous-models");
    }

    @Test
    public void testRootMeanSquareDeviationOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .command(CommandEnum.Measure.RMSD);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/expected/continuous-models");
    }

    @Test
    public void testRootMeanSquareDeviationOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .command(CommandEnum.Measure.RMSD);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/expected/incontinuous-models");
    }

    @Test
    public void testRootMeanSquareDeviationOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .command(CommandEnum.Measure.RMSD);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ROOT-MEAN-SQUARE-DEVIATION/expected/incontinuous-models");
    }

    @Test
    public void testPvalueOfSingleContinousModelWithinTheContextOfReferenceStructure() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .command(CommandEnum.Measure.P_VALUE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/expected/continuous-models");
    }

    @Test
    public void testPvalueOfMultipleContinousModelsWithinTheContextOfReferenceStructure() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .command(CommandEnum.Measure.P_VALUE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/expected/continuous-models");
    }

    @Test
    public void testPvalueOfSingleIncontinousModelWithinTheContextOfReferenceStructure() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .command(CommandEnum.Measure.P_VALUE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/expected/incontinuous-models");
    }

    @Test
    public void testPvalueOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .command(CommandEnum.Measure.P_VALUE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/P-VALUE/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityWatsonCrickOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_WC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityWatsonCrickOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_WC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityWatsonCrickOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_WC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityWatsonCrickOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_WC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-WATSON-CRICK/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityNonWatsonCrickOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_NWC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityNonWatsonCrickOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_NWC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityNonWatsonCrickOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_NWC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityNonWatsonCrickOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_NWC);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityStackingOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_STACKING);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityStackingOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_STACKING);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityStackingOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_STACKING);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityStackingOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_STACKING);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-STACKING/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityAllOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_ALL);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityAllOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_ALL);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/expected/continuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityAllOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_ALL);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/expected/incontinuous-models");
    }

    @Test
    public void testInteractionNetworkFidelityAllOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF_ALL);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/INTERACTION-NETWORK-FIDELITY-ALL/expected/incontinuous-models");
    }

    @Test
    public void testAllInteractionNetworkFidelityScoresAtOnceOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/expected/continuous-models");
    }

    @Test
    public void testAllInteractionNetworkFidelityScoresAtOnceOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/expected/continuous-models");
    }

    @Test
    public void testAllInteractionNetworkFidelityScoresAtOnceOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/expected/incontinuous-models");
    }

    @Test
    public void testAllInteractionNetworkFidelityScoresAtOnceOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.INF);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE/expected/incontinuous-models");
    }

    @Test
    public void testDeformationIndexOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.DI);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/expected/continuous-models");
    }

    @Test
    public void testDeformationIndexOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.DI);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/expected/continuous-models");
    }

    @Test
    public void testDeformationIndexOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.DI);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/expected/incontinuous-models");
    }

    @Test
    public void testDeformationIndexOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.DI);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/DEFORMATION-INDEX/expected/incontinuous-models");
    }

    @Test
    public void testAllScoresAtOnceOfSingleContinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/continuous-models/8_0_solution_4L81_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/continuous-models/models/8_Bujnicki_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.ALL_SCORES_AT_ONCE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/expected/continuous-models");
    }

    @Test
    public void testAllScoresAtOnceOfMultipleContinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/continuous-models/8_0_solution_4L81_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/continuous-models/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("8_Bujnicki_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.ALL_SCORES_AT_ONCE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/expected/continuous-models");
    }

    @Test
    public void testAllScoresAtOnceOfSingleIncontinousModelWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/incontinuous-models/models/14_ChenPostExp_1_rpr.pdb")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.ALL_SCORES_AT_ONCE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/expected/incontinuous-models");
    }

    @Test
    public void testAllScoresAtOnceOfMultipleIncontinousModelsWithinTheContextOfReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/incontinuous-models/14_5ddp_bound_solution_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/incontinuous-models/models")
                .targetModelsAlignment(
                        "14_5ddp_bound_solution_rpr.pdb:A_1,31|A_33,29;14_ChenPostExp_1_rpr.pdb:U_1,31|U_33,29")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("14_ChenPostExp_1_rpr.xml")
                .basePairsIdentificationTool(CommandEnum.BpsIdentificationTool.MC_ANNOTATE)
                .command(CommandEnum.Measure.ALL_SCORES_AT_ONCE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SCORES/ALL-SCORES-AT-ONCE/expected/incontinuous-models");
    }

    @Test
    public void testWholeSequenceOfSingleContinousModelWithSequenceIssues() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/SEQUENCE/1_solution_0_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/SEQUENCE/models/1_santalucia_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("1_santalucia_1_rpr.xml")
                .command(CommandEnum.Sequence.SEQUENCE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/SEQUENCE/expected");
    }

    @Test
    public void testWholeSequenceOfMultipleContinousModelsWithSequenceIssues() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/SEQUENCE/1_solution_0_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/SEQUENCE/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("1_santalucia_1_rpr.xml")
                .command(CommandEnum.Sequence.SEQUENCE);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/SEQUENCE/expected");
    }

    @Test
    public void testFragmentsOfSequenceOfSingleContinousModelWithSequenceIssues() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/FRAGMENT/1_solution_0_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/FRAGMENT/models/1_santalucia_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true)
                .outputFilePath("1_santalucia_1_rpr.xml")
                .targetModelsAlignment(
                        "1_solution_0_rpr.pdb:A_1,23|B_1,23;1_santalucia_1_rpr.pdb:A_1,23|B_1,23")
                .command(CommandEnum.Sequence.FRAGMENT);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/FRAGMENT/expected");
    }

    @Test
    public void testFragmentsOfSequenceOfMultipleContinousModelsWithSequenceIssues() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/FRAGMENT/1_solution_0_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/FRAGMENT/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true)
                .outputFilePath("1_santalucia_1_rpr.xml")
                .targetModelsAlignment(
                        "1_solution_0_rpr.pdb:A_1,23|B_1,23;1_santalucia_1_rpr.pdb:A_1,23|B_1,23")
                .command(CommandEnum.Sequence.FRAGMENT);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/SEQUENCE/FRAGMENT/expected");
    }

    @Test
    public void testOriginal3dOfSingleIncontinousModelSuperimposedOverReferenceStructure() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/ORIGINAL-3D/1_solution_0_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/ORIGINAL-3D/models/1_santalucia_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true)
                .outputFilePath("archive.zip")
                .targetModelsAlignment(
                        "1_solution_0_rpr.pdb:A_1,14|A_16,2|A_19,5|B_1,14|B_16,2|B_19,5;1_santalucia_1_rpr.pdb:A_1,14|A_16,2|A_19,5|B_1,14|B_16,2|B_19,5")
                .command(CommandEnum.Structure3d.ORIGINAL_3D);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/ORIGINAL-3D/expected");
    }

    @Test
    public void testOriginal3dOfMultipleIncontinousModelsSuperimposedOverReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/ORIGINAL-3D/1_solution_0_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/ORIGINAL-3D/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true)
                .outputFilePath("archive.zip")
                .targetModelsAlignment(
                        "1_solution_0_rpr.pdb:A_1,14|A_16,2|A_19,5|B_1,14|B_16,2|B_19,5;1_santalucia_1_rpr.pdb:A_1,14|A_16,2|A_19,5|B_1,14|B_16,2|B_19,5")
                .command(CommandEnum.Structure3d.ORIGINAL_3D);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/ORIGINAL-3D/expected");
    }

    @Test
    public void testRenumerated3dOfSingleContinousModelSuperimposedOverReferenceStructure() throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/RENUMERATED-3D/4_0_solution_3V7E_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/RENUMERATED-3D/models/4_adamiak_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .targetModelsAlignment("4_0_solution_3V7E_rpr.pdb:C_1,126;4_adamiak_1_rpr.pdb:A_1,126")
                .command(CommandEnum.Structure3d.RENUMERATED_3D);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/RENUMERATED-3D/expected");
    }

    @Test
    public void testRenumerated3dOfMultipleContinousModelSuperimposedOverReferenceStructure()
            throws Exception {
        final ComparisonInputModelImpl.Builder comparisonInputModelBuilder = new ComparisonInputModelImpl.Builder();
        comparisonInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/RENUMERATED-3D/4_0_solution_3V7E_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/RENUMERATED-3D/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .targetModelsAlignment("4_0_solution_3V7E_rpr.pdb:C_1,126;4_adamiak_1_rpr.pdb:A_1,126")
                .command(CommandEnum.Structure3d.RENUMERATED_3D);
        compare(comparisonInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/STRUCTURE-3D/RENUMERATED-3D/expected");
    }

    @Test
    public void testBasicDeformationProfileBasedOnSingleModel() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/target.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/models/model.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("model.svg")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/basic");
    }

    @Test
    public void testBasicDeformationProfileBasedOnMultipleModels() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/target.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/basic");
    }

    @Test
    public void testDataOnlyFromDeformationProfileBasedOnSingleModel() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .dataGenerated(true)
                .imageGenerated(false)
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/target.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/models/model.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/data");
    }

    @Test
    public void testDataOnlyFromDeformationProfileBasedOnMultipleModels() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .dataGenerated(true)
                .imageGenerated(false)
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/target.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/data");
    }

    @Test
    public void testCustomDeformationProfileBasedOnSingleModel() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .helices("H1,2,3,7,2;H2,12,2,16,2")
                .loops("L1,0,2;L2,5,2;L3,9,3;L4,14,2;L5,18,2")
                .draw("H1,H2,L1,L2,L3,L4,L5,H1xH2:I x II,L1xL2,H1xL2")
                .dataGenerated(true)
                .imageGenerated(true)
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/target.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/models/model.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/custom");
    }

    @Test
    public void testCustomDeformationProfileBasedOnMultipleModels() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .helices("H1,2,3,7,2;H2,12,2,16,2")
                .loops("L1,0,2;L2,5,2;L3,9,3;L4,14,2;L5,18,2")
                .draw("H1,H2,L1,L2,L3,L4,L5,H1xH2:I x II,L1xL2,H1xL2")
                .dataGenerated(true)
                .imageGenerated(true)
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/target.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/common/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/custom");
    }

    @Test
    public void testBasicDeformationProfileDrivenByAlignmentBasedOnSingleModel() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/alignment/4_0_solution_3V7E_rpr.pdb")
                .singlePdbModelFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/alignment/models/4_adamiak_1_rpr.pdb")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("4_adamiak_1_rpr.svg")
                .targetModelsAlignment("4_0_solution_3V7E_rpr.pdb:C_1,126;4_adamiak_1_rpr.pdb:A_1,126")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/basic-by-alignment");
    }

    @Test
    public void testBasicDeformationProfileDrivenByAlignmentBasedOnMultipleModels() throws Exception {
        final DeformationProfileInputModelImpl.Builder deformationProfileInputModelBuilder = new DeformationProfileInputModelImpl.Builder();
        deformationProfileInputModelBuilder
                .referenceStructurePdbFilePath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/alignment/4_0_solution_3V7E_rpr.pdb")
                .multiplePdbModelsDirectoryPath(
                        "ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/alignment/models")
                .considerAtomsSupportedByRNAPuzzlesOnly(true).outputFilePath("archive.zip")
                .targetModelsAlignment("4_0_solution_3V7E_rpr.pdb:C_1,126;4_adamiak_1_rpr.pdb:A_1,126")
                .command(CommandEnum.DEFORMATION_PROFILE);
        compare(deformationProfileInputModelBuilder,
                "./ANALYSIS-WITHIN-CONTEXT-OF-REFERENCE-STRUCTURE/DEFORMATION-PROFILE/expected/basic-by-alignment");
    }

    private void compare(final ComparisonInputModelImpl.Builder comparisonInputModelBuilder,
            final String expectedDirectory) throws Exception {
        final Class<?> clazz = this.getClass();
        final File referenceStructureFile = getFile(clazz, "./", "/",
                comparisonInputModelBuilder.getReferenceStructurePdbFilePath());
        comparisonInputModelBuilder.referenceStructurePdbFilePath(referenceStructureFile.getCanonicalPath());
        analyse(comparisonInputModelBuilder, expectedDirectory);
    }

    private void analyse(final AnalysisInputModelImpl.Builder analysisInputModelBuilder,
            final String expectedDirectory) throws Exception {
        final Class<?> clazz = this.getClass();
        File inputFile;
        if (StringUtils.isNotBlank(analysisInputModelBuilder.getSinglePdbModelFilePath())) {
            inputFile = getFile(clazz, "./", "/", analysisInputModelBuilder.getSinglePdbModelFilePath());
            analysisInputModelBuilder.singlePdbModelFilePath(inputFile.getCanonicalPath());
        } else if (StringUtils.isNotBlank(analysisInputModelBuilder.getMultiplePdbModelsDirectoryPath())) {
            inputFile = getFile(clazz, "./", "/",
                    analysisInputModelBuilder.getMultiplePdbModelsDirectoryPath());
            analysisInputModelBuilder.multiplePdbModelsDirectoryPath(inputFile.getCanonicalPath());
        }
        final String outputFileName = analysisInputModelBuilder.getOutputFilePath();
        final File outputFile = temporaryFolder.newFile(outputFileName);
        analysisInputModelBuilder.outputFilePath(outputFile.getCanonicalPath());
        final File expectedFile = getFile(clazz, expectedDirectory, "/", outputFileName);
        final AnalysisInputModel analysisInputModel = analysisInputModelBuilder.build();
        final String[] args = analysisInputModel.getArgs();
        assessment.perform(args);
        if (ArchiverFactory.isArchive(expectedFile)) {
            assertArchiveEquals(outputFile, expectedFile);
        } else {
            FileAssert.assertEquals(expectedFile, outputFile);
        }
    }

    private void assertArchiveEquals(final File outputFile, final File expectedFile) throws IOException {
        final Archiver archiver = ArchiverFactory.getArchiver(expectedFile);
        final File expectedFolder = temporaryFolder.newFolder();
        archiver.extract(expectedFile, expectedFolder);
        final File outputFolder = temporaryFolder.newFolder();
        archiver.extract(outputFile, outputFolder);
        final File[] expectedFolderFiles = expectedFolder.listFiles();
        final File[] outputFolderFiles = outputFolder.listFiles();
        final int expectedFolderFilesNo = ArrayUtils.getLength(expectedFolderFiles);
        final int outputFolderFilesNo = ArrayUtils.getLength(outputFolderFiles);
        Assert.assertEquals(expectedFolderFilesNo, outputFolderFilesNo);
        final Map<String, File> conversionResultFiles = Maps.newHashMap();
        for (File conversionResultFile : outputFolderFiles) {
            conversionResultFiles.put(conversionResultFile.getName(), conversionResultFile);
        }
        for (File file : expectedFolder.listFiles()) {
            FileAssert.assertEquals(file, conversionResultFiles.get(file.getName()));
        }
    }

    private static File getFile(final Class<?> clazz, final String... options) {
        final StringBuilder path = new StringBuilder();
        final int optionsCount = ArrayUtils.getLength(options);
        for (int optionIndex = 0; optionIndex < optionsCount; optionIndex++) {
            path.append(options[optionIndex]);
            if (optionIndex < optionsCount - 1) {
                path.append("/");
            }
        }
        return FileUtils.toFile(clazz.getResource(path.toString()));
    }

}
