package edu.put.ma.templates;

import java.util.Arrays;
import java.util.List;

import lombok.Getter;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.Lists;

public class ResidueTemplateImpl implements ResidueTemplate {

    @Getter
    private String oneResidueName = null;

    @Getter
    private String threeResidueName = null;

    private List<AtomTemplate> atomTemplates = null;

    private List<String> alterResNames = null;

    public ResidueTemplateImpl(final String residueNames, final String atomNames,
            final String templateAtomName) {
        alterResNames = Arrays.asList(StringUtils.split(residueNames, ','));
        if (CollectionUtils.size(alterResNames) >= 2) {
            oneResidueName = alterResNames.get(0);
            threeResidueName = alterResNames.get(1);
        }
        atomTemplates = Lists.newArrayList(Arrays.asList(new AtomTemplate[] { new AtomTemplateImpl(atomNames,
                templateAtomName) }));
    }

    @Override
    public final void addAtomTemplate(final String atomNames, final String templateAtomName) {
        atomTemplates.add(new AtomTemplateImpl(atomNames, templateAtomName));
    }

    @Override
    public final boolean isEqualResidueName(final String presName) {
        final String resName = presName.trim().toUpperCase();
        if ((resName.equals(threeResidueName)) || (resName.equals(oneResidueName))
                || (resName.endsWith(oneResidueName))) {
            return true;
        }
        return false;
    }

    @Override
    public final boolean isProperAtomName(final String atomName) {
        for (int i = 0; i < atomTemplates.size(); i++) {
            if (atomTemplates.get(i).isProperAtomName(atomName)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public final String getTemplateAtomName(final String atomName) {
        return atomTemplates.get(getProperAtomNameIndex(atomName)).getTemplateAtomName();
    }

    @Override
    public final String getOriginalAtomName(final String atomName) {
        return atomTemplates.get(getProperAtomNameIndex(atomName)).getOriginalAtomName();
    }

    @Override
    public boolean isAppropriateResName(final String resName) {
        return alterResNames.contains(resName);
    }

    private int getProperAtomNameIndex(final String atomName) {
        for (int i = 0; i < atomTemplates.size(); i++) {
            if (atomTemplates.get(i).isProperAtomName(atomName)) {
                return i;
            }
        }
        throw new IllegalArgumentException(String.format("No atom %s found!", atomName));
    }
}
