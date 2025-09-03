package edu.put.ma.model.input;

import java.util.List;
import java.util.regex.Pattern;

import lombok.Getter;

import org.apache.commons.cli.Options;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;

import edu.put.ma.archiver.ArchiverType;
import edu.put.ma.model.services.DeformationProfile;
import edu.put.ma.model.services.DpParameter;
import edu.put.ma.model.services.DrawComponent;
import edu.put.ma.model.services.Helix;
import edu.put.ma.model.services.Loop;
import edu.put.ma.model.services.StructuresSet;
import static edu.put.ma.Assessment.SUPPORTED_IMAGE_EXTENSION;

public class DeformationProfileInputModelImpl extends ComparisonInputModelImpl {

    private static final String HELIX_LABEL = "helix";

    private static final String LOOP_LABEL = "loop";

    private static final Logger LOGGER = LoggerFactory.getLogger(DeformationProfileInputModelImpl.class);

    private String helices;

    private String loops;

    private String draw;

    private boolean dataGenerated = false;

    private boolean imageGenerated = true;

    public DeformationProfileInputModelImpl(final String[] args, final String artifactId) {
        super(args, artifactId);
        secureInitState(artifactId);
    }

    protected DeformationProfileInputModelImpl(final Builder deformationProfileInputModelBuilder) {
        super(deformationProfileInputModelBuilder);
        this.helices = deformationProfileInputModelBuilder.helices;
        this.loops = deformationProfileInputModelBuilder.loops;
        this.draw = deformationProfileInputModelBuilder.draw;
        this.dataGenerated = deformationProfileInputModelBuilder.dataGenerated;
        this.imageGenerated = deformationProfileInputModelBuilder.imageGenerated;
        this.initOptionsMapping();
    }

    @Override
    public boolean isInputInitializedProperly() {
        return super.isInputInitializedProperly() && isValid();
    }

    @Override
    public Options constructSpecificOptions() {
        final Options options = super.constructSpecificOptions();
        options.addOption("h", "helices", true, "(optional) helices");
        options.addOption("l", "loops", true, "(optional) loops");
        options.addOption("w", "draw", true, "(optional) draw");
        options.addOption("m", "generate-matrix-data", false, "(optional) generate matrix data");
        options.addOption("g", "generate-svg-image", false, "(optional) generate svg image");
        return options;
    }

    @Override
    public boolean isValid() {
        return isSecondaryStructureStringValid(helices, HELIX_LABEL)
                && isSecondaryStructureStringValid(loops, LOOP_LABEL) && isDrawStringValid(draw);
    }

    @Override
    public StructuresSet constructStructuresSet() {
        final StructuresSet structuresSet = super.constructStructuresSet();
        structuresSet.setDeformationProfile(this.getDeformationProfile());
        return structuresSet;
    }

    @Override
    public String getInputModelString() {
        final StringBuilder inputModelStringBuilder = new StringBuilder(super.getInputModelString());
        edu.put.ma.utils.StringUtils.extendStringBuilder(inputModelStringBuilder, helices, "Helices");
        edu.put.ma.utils.StringUtils.extendStringBuilder(inputModelStringBuilder, loops, "Loops");
        edu.put.ma.utils.StringUtils.extendStringBuilder(inputModelStringBuilder, draw, "Draw");
        inputModelStringBuilder.append("\nGenerate matrix data: ").append(dataGenerated ? "Y" : "N");
        inputModelStringBuilder.append("\nGenerate svg image: ").append(imageGenerated ? "Y" : "N");
        return inputModelStringBuilder.toString();
    }

    @Override
    protected void initOptionsMapping() {
        super.initOptionsMapping();
        optionsMapping.putAll(new ImmutableMap.Builder<String, String>().put("helices", "-h")
                .put("loops", "-l").put("draw", "-w").put("dataGenerated", "-m").put("imageGenerated", "-g")
                .build());
    }

    @Override
    protected void initState() {
        super.initState();
        setHelices();
        setLoops();
        setDraw();
        setDataGenerated();
        setImageGenerated();
        imageGenerated = ((!dataGenerated) && (!imageGenerated)) ? true : imageGenerated;
        refineOutputFilePath();
    }

    @Getter
    public static class Builder extends ComparisonInputModelImpl.Builder {

        private String helices;

        private String loops;

        private String draw;

        private boolean dataGenerated = false;

        private boolean imageGenerated = true;

        public Builder helices(final String helices) {
            this.helices = edu.put.ma.utils.StringUtils.deleteAllSpaces(helices);
            return this;
        }

        public Builder loops(final String loops) {
            this.loops = edu.put.ma.utils.StringUtils.deleteAllSpaces(loops);
            return this;
        }

        public Builder draw(final String draw) {
            this.draw = StringUtils.trim(draw);
            return this;
        }

        public Builder dataGenerated(final boolean dataGenerated) {
            this.dataGenerated = dataGenerated;
            return this;
        }

        public Builder imageGenerated(final boolean imageGenerated) {
            this.imageGenerated = imageGenerated;
            return this;
        }

        public AnalysisInputModel build() {
            return new DeformationProfileInputModelImpl(this);
        }
    }

    private void setHelices() {
        helices = getOptionString("h");
    }

    private void setLoops() {
        loops = getOptionString("l");
    }

    private void setDraw() {
        draw = getOptionString("w");
    }

    private void setDataGenerated() {
        dataGenerated = isOptionSet("m");
    }

    private void setImageGenerated() {
        imageGenerated = isOptionSet("g");
    }

    private DeformationProfile getDeformationProfile() {
        final DeformationProfile profile = new DeformationProfile();
        profile.addAllLoops(getLoops(loops));
        profile.addAllHelices(getHelices(helices));
        profile.addAllDraws(getDraws(draw, profile.getHelicesList(), profile.getLoopsList()));
        profile.setIsDataGenerated(dataGenerated);
        profile.setIsImageGenerated(imageGenerated);
        return profile;
    }

    private void refineOutputFilePath() {
        final String outputFileExtension = StringUtils.lowerCase(FilenameUtils
                .getExtension(getOutputFilePath()));
        if (isMultipleModelsDirectoryAtInput()
                && (!ArchiverType.ZIP.getExtension().equals(outputFileExtension))) {
            refineOutputFilePath(outputFileExtension, ArchiverType.ZIP.getExtension());
        } else if (dataGenerated && (!"zip".equals(outputFileExtension))) {
            refineOutputFilePath(outputFileExtension, ArchiverType.ZIP.getExtension());
        } else {
            final boolean isImage = !dataGenerated && imageGenerated;
            if (isSingleModelAtInput() && isImage && (!SUPPORTED_IMAGE_EXTENSION.equals(outputFileExtension))) {
                refineOutputFilePath(outputFileExtension, SUPPORTED_IMAGE_EXTENSION);
            }
        }
    }

    private void refineOutputFilePath(final String outputFileExtension, final String expectedExtension) {
        if (StringUtils.isBlank(outputFileExtension)) {
            setOutputFilePath(new StringBuilder(getOutputFilePath()).append(".").append(expectedExtension)
                    .toString());
        } else {
            setOutputFilePath(getOutputFilePath().replaceAll(
                    new StringBuilder("\\.").append(outputFileExtension).toString(),
                    new StringBuilder(".").append(expectedExtension).toString()));
        }
    }

    private static List<Loop> getLoops(final String loopsString) {
        final List<Loop> loops = Lists.newArrayList();
        setSecondaryStructures(loopsString, loops, false, LOOP_LABEL);
        return loops;
    }

    private static List<Helix> getHelices(final String helicesString) {
        final List<Helix> helices = Lists.newArrayList();
        setSecondaryStructures(helicesString, helices, true, HELIX_LABEL);
        return helices;
    }

    private static List<DrawComponent> getDraws(final String drawString, final List<Helix> helices,
            final List<Loop> loops) {
        final List<DrawComponent> result = Lists.newArrayList();
        if ((StringUtils.isNotBlank(drawString)) && (StringUtils.indexOf(drawString, ',') >= 0)) {
            final String[] elements = StringUtils.split(drawString, ',');
            for (String element : elements) {
                if (StringUtils.indexOf(element, ':') >= 0) {
                    final String[] compound = StringUtils.split(element, ':');
                    parseCompound(helices, loops, result, compound);
                } else {
                    parseCompound(helices, loops, result, new String[] { element });
                }
            }
        }
        return result;
    }

    private static <T extends DpParameter> void setSecondaryStructures(final String secondaryStructureString,
            final List<T> secondaryStructures, final boolean isHelix, final String type) {
        if (StringUtils.isNotBlank(secondaryStructureString)) {
            if (StringUtils.indexOf(secondaryStructureString, ';') >= 0) {
                final String[] currentHelices = StringUtils.split(secondaryStructureString, ';');
                analyseFragments(secondaryStructures, isHelix, type, currentHelices);
            } else {
                if (isFragmentCorrect(1, secondaryStructureString, type)) {
                    addSecondaryStructure(secondaryStructureString, secondaryStructures, isHelix, type);
                }
            }
        }
    }

    private static <T extends DpParameter> void analyseFragments(final List<T> secondaryStructures,
            final boolean isHelix, final String type, final String[] currentHelices) {
        int fragmentNo = 1;
        for (String fragment : currentHelices) {
            if (isFragmentCorrect(fragmentNo, fragment, type)) {
                addSecondaryStructure(fragment, secondaryStructures, isHelix, type);
            }
            fragmentNo++;
        }
    }

    @SuppressWarnings("unchecked")
    private static <T extends DpParameter> void addSecondaryStructure(final String fragmentString,
            final List<T> fragments, final boolean isHelix, final String type) {
        final String[] elements = StringUtils.split(fragmentString, ',');
        if (elements.length == ((isHelix) ? 5 : 3)) {
            T fragment;
            if (isHelix) {
                fragment = (T) new Helix(elements);
            } else {
                fragment = (T) new Loop(elements);
            }
            if (!fragments.contains(fragment)) {
                fragments.add(fragment);
            } else {
                final String fString = (isHelix) ? ((Helix) fragment).toString() : ((Loop) fragment)
                        .toString();
                LOGGER.info(String
                        .format("%s will be omitted because there was other %s defined previously with the same name (%s)!",
                                fString, type, elements[0]));
            }
        } else {
            LOGGER.info(String
                    .format("%s will be omitted because it is not compliant with the following scheme: name + ',' + resIdx5 + ',' + distance + ',' + resIdx3 + ',' + distance",
                            fragmentString));
        }
    }

    private static <T extends DpParameter> boolean secondaryStructureExistsWithName(
            final String elementWithoutSpaces, final List<T> secondaryStructures) {
        if ((StringUtils.isNotBlank(elementWithoutSpaces)) && (secondaryStructures != null)) {
            for (T secondaryStructure : secondaryStructures) {
                if (StringUtils.equalsIgnoreCase(elementWithoutSpaces, secondaryStructure.getName())) {
                    return true;
                }
            }
        }
        return false;
    }

    private static void parseCompound(final List<Helix> helices, final List<Loop> loops,
            final List<DrawComponent> result, final String[] compound) {
        final String elementWithoutSpaces = edu.put.ma.utils.StringUtils.deleteAllSpaces(compound[0]);
        if (StringUtils.indexOf(elementWithoutSpaces, 'x') >= 0) {
            final String[] names = StringUtils.split(elementWithoutSpaces, 'x');
            if (names.length == PAIR_SIZE) {
                if (ArrayUtils.getLength(compound) == PAIR_SIZE) {
                    addDrawComponent(helices, loops, result, names, elementWithoutSpaces, compound[1]);
                } else {
                    addDrawComponent(helices, loops, result, names, elementWithoutSpaces);
                }
            } else {
                LOGGER.info(String.format(
                        "%s will be omitted because the number of considered components is various than 2.",
                        elementWithoutSpaces));
            }
        } else {
            addDrawComponent(helices, loops, result, new String[] { elementWithoutSpaces },
                    elementWithoutSpaces);
        }
    }

    private static boolean secondaryStructureExistsWithName(final String[] names, final List<Helix> helices,
            final List<Loop> loops) {
        boolean result = true;
        for (String name : names) {
            if (!result) {
                break;
            }
            result = result
                    && ((secondaryStructureExistsWithName(name, helices)) || (secondaryStructureExistsWithName(
                            name, loops)));
        }
        return result;
    }

    private static void addDrawComponent(final List<Helix> helices, final List<Loop> loops,
            final List<DrawComponent> result, final String[] names, final String element) {
        addDrawComponent(helices, loops, result, names, element, null);
    }

    private static void addDrawComponent(final List<Helix> helices, final List<Loop> loops,
            final List<DrawComponent> result, final String[] names, final String element,
            final String description) {
        final String name = getName(names);
        if (secondaryStructureExistsWithName(names, helices, loops)) {
            final DrawComponent newComponent = (StringUtils.isBlank(description)) ? new DrawComponent(element)
                    : new DrawComponent(element, description);
            if (!result.contains(newComponent)) {
                result.add(newComponent);
            } else {
                LOGGER.info(String
                        .format("%s will be omitted because there was other secondary structure defined earlier with the same name (%s).",
                                newComponent.toString(), element));
            }
        } else {
            LOGGER.info(String.format("There is no either helix or loop with name %s.", name));
        }
    }

    private static String getName(final String[] names) {
        final StringBuilder nameSb = new StringBuilder();
        for (String name : names) {
            if (nameSb.length() > 0) {
                nameSb.append(", ");
            }
            nameSb.append(name);
        }
        return nameSb.toString();
    }

    private static boolean isSecondaryStructureStringValid(final String secondaryStructureString,
            final String type) {
        boolean isOk = true;
        if (StringUtils.isNotBlank(secondaryStructureString)) {
            if (StringUtils.indexOf(secondaryStructureString, ';') >= 0) {
                final String[] currentFragments = StringUtils.split(secondaryStructureString, ';');
                int fragmentNo = 1;
                for (String fragment : currentFragments) {
                    isOk = isFragmentCorrect(fragmentNo++, fragment, type);
                }
            } else {
                isOk = isFragmentCorrect(1, secondaryStructureString, type);
            }
        }
        return isOk;
    }

    private static boolean isFragmentCorrect(final int fragmentNo, final String fragment, final String type) {
        boolean isOk = true;
        final String[] elements = StringUtils.split(fragment, ',');
        for (int i = 1; i < elements.length; i++) {
            final String val = StringUtils.trim(elements[i]);
            if (i % 2 == 1) {
                if (!Pattern.matches("[-]?[0-9]+", val)) {
                    LOGGER.info(String.format("ResNo (%s) of %s (%d) should be integer number e.g., 5.",
                            type, val, fragmentNo));
                    isOk = false;
                }
            } else {
                if (!Pattern.matches("[0-9]+", val)) {
                    LOGGER.info(String.format("Distance (%s) of %s (%d) should be natural number e.g., 5.",
                            type, val, fragmentNo));
                    isOk = false;
                }
            }
        }
        return isOk;
    }

    private static boolean isDrawStringValid(final String drawString) {
        boolean isOk = true;
        if (StringUtils.isNotBlank(drawString)) {
            List<String> singleElements = Lists.newArrayList();
            List<String> singleElementNames = Lists.newArrayList();
            List<String> compounds = Lists.newArrayList();
            List<String> compoundNames = Lists.newArrayList();
            isOk = parseDrawString(drawString, singleElements, singleElementNames, compounds, compoundNames);
            int compoundIdx = 0;
            for (String compound : compounds) {
                isOk = verifyCompoundOfDrawString(singleElements, singleElementNames,
                        compoundNames.get(compoundIdx++), compound);
            }
        }
        return isOk;
    }

    private static boolean verifyCompoundOfDrawString(final List<String> singleElements,
            final List<String> singleElementNames, final String compoundName, final String compound) {
        boolean isOk = true;
        final String[] elements = StringUtils.split(compound, 'x');
        for (int i = 0; i < elements.length; i++) {
            if (!singleElements.contains(elements[i])) {
                LOGGER.info(String.format(
                        "Secondary structure component (%s) included in the compound (%s) is not found!",
                        elements[i], compound));
                isOk = false;
            }
            final int fIdx = singleElements.indexOf(elements[i]);
            if ((!StringUtils.equals("", singleElementNames.get(fIdx)))
                    && (StringUtils.indexOf(compoundName, singleElementNames.get(fIdx)) < 0)) {
                LOGGER.info("Structure component (" + compound + ") description (" + compoundName
                        + ") does not contain description of " + elements[i] + "("
                        + singleElementNames.get(fIdx) + ")!");
            }
        }
        return isOk;
    }

    private static boolean parseDrawString(final String drawString, final List<String> singleElements,
            final List<String> singleElementNames, final List<String> compounds,
            final List<String> compoundNames) {
        boolean isOk = true;
        if (StringUtils.indexOf(drawString, ',') >= 0) {
            final String[] elements = StringUtils.split(drawString, ',');
            for (String element : elements) {
                if (StringUtils.indexOf(element, ':') >= 0) {
                    final String[] compound = StringUtils.split(element, ':');
                    isOk = isOk
                            && analyseCompound(singleElements, singleElementNames, compounds, compoundNames,
                                    compound);
                } else {
                    analyseSingleElementOfDrawString(singleElements, singleElementNames, compounds,
                            compoundNames, element, "");
                }
            }
        }
        return isOk;
    }

    private static boolean analyseCompound(final List<String> singleElements,
            final List<String> singleElementNames, final List<String> compounds,
            final List<String> compoundNames, final String[] compound) {
        if (ArrayUtils.getLength(compound) == PAIR_SIZE) {
            analyseSingleElementOfDrawString(singleElements, singleElementNames, compounds, compoundNames,
                    compound[0], compound[1]);
        } else {
            LOGGER.info("Secondary structure component can be described by a name or name:description e.g., H1 or H1:I.");
            return false;
        }
        return true;
    }

    private static void analyseSingleElementOfDrawString(final List<String> singleElements,
            final List<String> singleElementNames, final List<String> compounds,
            final List<String> compoundNames, final String element, final String name) {
        final String elementWithoutSpaces = edu.put.ma.utils.StringUtils.deleteAllSpaces(element);
        if (StringUtils.indexOf(elementWithoutSpaces, 'x') >= 0) {
            addComponent(compounds, compoundNames, name, elementWithoutSpaces);
        } else {
            addComponent(singleElements, singleElementNames, name, elementWithoutSpaces);
        }
    }

    private static void addComponent(final List<String> components, final List<String> componentNames,
            final String name, final String value) {
        if (!components.contains(value)) {
            components.add(value);
            componentNames.add(name);
        }
    }
}
