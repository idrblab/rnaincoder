package edu.put.ma.templates;

import java.util.List;

import lombok.Getter;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.Lists;

public class TemplatesImpl implements Templates {

    private static final int FIELDS_NO = 3;

    private List<ResidueTemplate> resTemplates = null;

    @Getter
    private String templatesString = null;

    @Override
    public final void load(final String dictionaryString) {
        if (this.resTemplates == null) {
            this.resTemplates = Lists.newArrayList();
        } else {
            this.resTemplates.clear();
        }
        this.templatesString = dictionaryString;
        ResidueTemplate resTemplate = null;
        for (String line : StringUtils.split(dictionaryString, '\n')) {
            if (StringUtils.indexOf(line, ';') >= 0) {
                final String[] elems = StringUtils.split(line, ';');
                if (elems.length == FIELDS_NO) {
                    resTemplate = constructResidueTemplate(resTemplate, elems);
                }
            }
        }
        if (resTemplate != null) {
            this.resTemplates.add(resTemplate);
        }
    }

    @Override
    public final String getThreeLetterResidueNameByStructuralName(final String resName) {
        final int index = this.getResidueIndexByName(resName);
        return (index >= 0) ? this.resTemplates.get(index).getThreeResidueName() : null;
    }

    @Override
    public final String getOneLetterResidueNameByStructuralName(final String resName) {
        final int index = this.getResidueIndexByName(resName);
        return (index >= 0) ? this.resTemplates.get(index).getOneResidueName() : null;
    }

    @Override
    public final String getTemplateAtomName(final String resName, final String atomName) {
        final int index = this.getResidueIndexByName(resName);
        return (index >= 0) ? this.resTemplates.get(index).getTemplateAtomName(atomName) : null;
    }

    @Override
    public final String getOriginalAtomName(final String resName, final String atomName) {
        final int index = this.getResidueIndexByName(resName);
        return (index >= 0) ? this.resTemplates.get(index).getOriginalAtomName(atomName) : null;
    }

    @Override
    public final boolean isProperAtomName(final String residueName, final String atomName) {
        final int index = this.getResidueIndexByName(residueName);
        return (index >= 0) ? this.resTemplates.get(index).isProperAtomName(atomName) : false;
    }

    @Override
    public final boolean isInitiated() {
        return this.resTemplates != null;
    }

    private ResidueTemplate constructResidueTemplate(final ResidueTemplate resTemplate, final String[] elems) {
        ResidueTemplate residueTemplate = resTemplate;
        if (residueTemplate == null) {
            residueTemplate = new ResidueTemplateImpl(elems[0], elems[1], elems[2]);
        } else {
            if (residueTemplate.isEqualResidueName(StringUtils.split(elems[0], ',')[1])) {
                residueTemplate.addAtomTemplate(elems[1], elems[2]);
            } else {
                this.resTemplates.add(residueTemplate);
                residueTemplate = new ResidueTemplateImpl(elems[0], elems[1], elems[2]);
            }
        }
        return residueTemplate;
    }

    private int getResidueIndexByName(final String pResName) {
        final String resName = pResName.trim().toUpperCase();
        if (resName.length() > 1) {
            return getResidueIndexBasedOnManyLetterName(resName);
        } else {
            return getResidueIndexBasedOnOneLetterName(resName);
        }
    }

    private int getResidueIndexBasedOnOneLetterName(final String resName) {
        for (int i = 0; i < this.resTemplates.size(); i++) {
            if (resName.equals(this.resTemplates.get(i).getOneResidueName())) {
                return i;
            }
        }
        return -1;
    }

    private int getResidueIndexBasedOnManyLetterName(final String resName) {
        for (int i = 0; i < this.resTemplates.size(); i++) {
            if (this.resTemplates.get(i).isAppropriateResName(resName)) {
                return i;
            }
        }
        int howMany = 0;
        int index = -1;
        for (int i = 0; i < this.resTemplates.size(); i++) {
            if (resName.endsWith(this.resTemplates.get(i).getOneResidueName())) {
                howMany++;
                index = i;
            }
        }
        if (howMany == 1) {
            return index;
        }
        return -1;
    }
}
