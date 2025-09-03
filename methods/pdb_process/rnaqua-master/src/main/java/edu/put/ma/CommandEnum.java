package edu.put.ma;

import java.util.EnumSet;
import java.util.List;

import lombok.AllArgsConstructor;
import lombok.Getter;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import edu.put.ma.utils.ArrayUtils;

public enum CommandEnum implements Command {

    PDB_VALIDATION(new Mapping("PDB-VALIDATION", "validate")), DEFORMATION_PROFILE(new Mapping(
            "DEFORMATION-PROFILE", "dp"));

    @Getter
    private Mapping mapping;

    CommandEnum(final Mapping mapping) {
        this.mapping = mapping;
    }

    @Override
    public String toString() {
        return mapping.name;
    }

    public enum BpsIdentificationTool implements Command {

        RNAVIEW(new Mapping("RNAVIEW", "RNAView")), MC_ANNOTATE(new Mapping("MC-ANNOTATE", "MC-Annotate"));

        @Getter
        private Mapping mapping;

        BpsIdentificationTool(final Mapping mapping) {
            this.mapping = mapping;
        }

        @Override
        public String toString() {
            return mapping.name;
        }
    }

    public enum Measure implements Command {

        CS(new Mapping("CLASH-SCORE", "clashscore")), RMSD(new Mapping("ROOT-MEAN-SQUARE-DEVIATION", "rmsd")), INF(
                new Mapping("ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE", "inf")), INF_WC(new Mapping(
                "INTERACTION-NETWORK-FIDELITY-WATSON-CRICK", "infwc")), INF_NWC(new Mapping(
                "INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK", "infnwc")), INF_STACKING(new Mapping(
                "INTERACTION-NETWORK-FIDELITY-STACKING", "infstacking")), INF_ALL(new Mapping(
                "INTERACTION-NETWORK-FIDELITY-ALL", "infall")), P_VALUE(new Mapping("P-VALUE", "pvalue")), DI(
                new Mapping("DEFORMATION-INDEX", "di")), ALL_SCORES_AT_ONCE(new Mapping("ALL-SCORES-AT-ONCE",
                "cs"));

        @Getter
        private Mapping mapping;

        Measure(final Mapping mapping) {
            this.mapping = mapping;
        }

        @Override
        public String toString() {
            return mapping.name;
        }

        public static final List<String> getReferenceStructureBasedMeasures() {
            final List<String> referenceStructureBasedMeasures = Lists.newArrayList(ArrayUtils
                    .getEnumNames(CommandEnum.Measure.class));
            referenceStructureBasedMeasures.remove(CommandEnum.Measure.CS.toString());
            return referenceStructureBasedMeasures;
        }

    }

    public enum Sequence implements Command {

        SEQUENCE(new Mapping("SEQUENCE", "sequence")), FRAGMENT(new Mapping("FRAGMENT", "fragmentsSequence"));

        @Getter
        private Mapping mapping;

        Sequence(final Mapping mapping) {
            this.mapping = mapping;
        }

        @Override
        public String toString() {
            return mapping.name;
        }
    }

    public enum Structure3d implements Command {

        ORIGINAL_3D(new Mapping("ORIGINAL-3D", "3do")), RENUMERATED_3D(new Mapping("RENUMERATED-3D", "3dr"));

        @Getter
        private Mapping mapping;

        Structure3d(final Mapping mapping) {
            this.mapping = mapping;
        }

        @Override
        public String toString() {
            return mapping.name;
        }
    }

    public static <T extends Enum<T>> T fromString(final String name, final Class<T> enumClass) {
        Preconditions.checkNotNull(name, "Command should be defined!");
        for (T command : EnumSet.allOf(enumClass)) {
            if ((command instanceof Command)
                    && (StringUtils.equalsIgnoreCase(name, ((Command) command).getMapping().name))) {
                return command;
            }
        }
        throw new IllegalArgumentException("No command found!");
    }

    @AllArgsConstructor
    public static class Mapping {
        @Getter
        private String name;

        @Getter
        private String serviceCode;
    }

    public static final List<String> getSequenceAndStructure3dServices() {
        final List<String> result = Lists.newArrayList(ArrayUtils.getEnumNames(CommandEnum.Sequence.class));
        result.addAll(ArrayUtils.getEnumNames(CommandEnum.Structure3d.class));
        return result;
    }
}
