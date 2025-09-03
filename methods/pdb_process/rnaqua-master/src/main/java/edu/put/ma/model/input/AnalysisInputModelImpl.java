package edu.put.ma.model.input;

import java.io.File;

import lombok.Getter;
import lombok.NoArgsConstructor;

import org.apache.commons.cli.Options;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.ImmutableMap;

import edu.put.ma.CommandEnum;
import edu.put.ma.model.services.StructuresSet;
import edu.put.ma.utils.ArrayUtils;
import edu.put.ma.utils.PdbReader;
import edu.put.ma.utils.PreconditionUtils;

public class AnalysisInputModelImpl extends CommonInputModelImpl implements AnalysisInputModel {

    private static final Logger LOGGER = LoggerFactory.getLogger(AnalysisInputModelImpl.class);

    private String singlePdbModelFilePath;

    private String multiplePdbModelsDirectoryPath;

    public AnalysisInputModelImpl(final String[] args, final String artifactId) {
        super(args);
        secureInitState(artifactId);
    }

    protected AnalysisInputModelImpl(final Builder analysisInputModelBuilder) {
        super(analysisInputModelBuilder);
        this.singlePdbModelFilePath = analysisInputModelBuilder.singlePdbModelFilePath;
        this.multiplePdbModelsDirectoryPath = analysisInputModelBuilder.multiplePdbModelsDirectoryPath;
        this.initOptionsMapping();
    }

    @Override
    public boolean isInputInitializedProperly() {
        return super.isInputInitializedProperly()
                && (isCommandLineHasOption("s") || isCommandLineHasOption("d")) && isValid();
    }

    @Override
    public Options constructSpecificOptions() {
        final Options options = new Options();
        options.addOption("s", "single-model-file-path", true, "single model PDB file path");
        options.addOption("d", "multiple-models-directory-path", true, "multiple PDB models directory path");
        return options;
    }

    @Override
    public boolean isValid() {
        return isAppropriateCommandUsed() && isInputModelsValid();
    }

    @Override
    public StructuresSet constructStructuresSet() {
        final StructuresSet structuresSet = super.constructStructuresSet();
        if (StringUtils.isNotBlank(singlePdbModelFilePath)) {
            PdbReader.readSingleStructure(singlePdbModelFilePath, structuresSet, false);
        } else {
            PdbReader.readMultiplePdbModelsDirectory(multiplePdbModelsDirectoryPath, structuresSet);
        }
        return structuresSet;
    }

    @Override
    public String getInputModelString() {
        final StringBuilder inputModelStringBuilder = new StringBuilder(super.getInputModelString());
        if ((StringUtils.isNotBlank(singlePdbModelFilePath))
                && (StringUtils.isBlank(multiplePdbModelsDirectoryPath))) {
            inputModelStringBuilder.append("Single model PDB file path: ").append(singlePdbModelFilePath);
        } else if ((StringUtils.isBlank(singlePdbModelFilePath))
                && (StringUtils.isNotBlank(multiplePdbModelsDirectoryPath))) {
            inputModelStringBuilder.append("Multiple PDB models directory path: ").append(
                    multiplePdbModelsDirectoryPath);
        }
        return inputModelStringBuilder.toString();
    }

    @Override
    protected void initOptionsMapping() {
        super.initOptionsMapping();
        optionsMapping.putAll(ImmutableMap.of("singlePdbModelFilePath", "-s",
                "multiplePdbModelsDirectoryPath", "-d"));
    }

    @Override
    protected void initState() {
        super.initState();
        setSinglePdbModelFilePath();
        setMultiplePdbModelsDirectoryPath();
    }

    @Getter
    @NoArgsConstructor
    public static class Builder extends CommonInputModelImpl.Builder {

        private String singlePdbModelFilePath;

        private String multiplePdbModelsDirectoryPath;

        public Builder singlePdbModelFilePath(final String singlePdbModelFilePath) {
            this.singlePdbModelFilePath = singlePdbModelFilePath;
            return this;
        }

        public Builder multiplePdbModelsDirectoryPath(final String multiplePdbModelsDirectoryPath) {
            this.multiplePdbModelsDirectoryPath = multiplePdbModelsDirectoryPath;
            return this;
        }

        public AnalysisInputModel build() {
            return new AnalysisInputModelImpl(this);
        }
    }

    protected boolean isMultipleModelsDirectoryAtInput() {
        return (StringUtils.isNotBlank(multiplePdbModelsDirectoryPath))
                && (StringUtils.isBlank(singlePdbModelFilePath));
    }

    protected boolean isSingleModelAtInput() {
        return (StringUtils.isBlank(multiplePdbModelsDirectoryPath))
                && (StringUtils.isNotBlank(singlePdbModelFilePath));
    }

    private boolean isAppropriateCommandUsed() {
        return CommandEnum.Measure.CS == getCommand() || CommandEnum.PDB_VALIDATION == getCommand()
                || ArrayUtils.getEnumNames(CommandEnum.Sequence.class).contains(getCommand().toString())
                || ArrayUtils.getEnumNames(CommandEnum.Structure3d.class).contains(getCommand().toString());
    }

    private boolean isInputModelsValid() {
        if ((StringUtils.isBlank(singlePdbModelFilePath))
                && (StringUtils.isBlank(multiplePdbModelsDirectoryPath))) {
            LOGGER.info("You should define single model PDB file path (using '-s' option) or multiple PDB models directory path (using '-d' option).");
            return false;
        } else if ((StringUtils.isNotBlank(singlePdbModelFilePath))
                && (StringUtils.isNotBlank(multiplePdbModelsDirectoryPath))) {
            LOGGER.info("You should choose one option only ('-s' or '-d').");
            return false;
        } else if (StringUtils.isNotBlank(singlePdbModelFilePath)) {
            final File singlePdbModelFile = FileUtils.getFile(singlePdbModelFilePath);
            PreconditionUtils
                    .checkIfFileExistsAndIsNotADirectory(singlePdbModelFile, "Single PDB model file");
            return StringUtils.endsWithIgnoreCase(singlePdbModelFilePath, ".pdb");
        } else if (StringUtils.isNotBlank(multiplePdbModelsDirectoryPath)) {
            final File multiplePdbModelsDirectory = FileUtils.getFile(multiplePdbModelsDirectoryPath);
            PreconditionUtils.checkIfDirectoryExists(multiplePdbModelsDirectory,
                    "Multiple PDB models directory");
        }
        return true;
    }

    private void setSinglePdbModelFilePath() {
        singlePdbModelFilePath = getOptionString("s");
    }

    private void setMultiplePdbModelsDirectoryPath() {
        multiplePdbModelsDirectoryPath = getOptionString("d");
    }
}
