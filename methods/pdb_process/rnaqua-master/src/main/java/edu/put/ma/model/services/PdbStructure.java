package edu.put.ma.model.services;

import java.util.List;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlSeeAlso;

@XmlRootElement(name = "structure")
@XmlSeeAlso({ Fragment.class })
public class PdbStructure extends Description {

    @XmlElementWrapper(name = "alignment")
    @XmlElement(name = "fragment")
    private List<Fragment> fragmentList = null;

    @XmlElementWrapper(name = "atoms")
    @XmlElement(name = "atom")
    private List<String> atomList = null;

    public final List<Fragment> getFragmentsList() {
        return fragmentList;
    }

    public final void setFragmentList(List<Fragment> pfragmentList) {
        this.fragmentList = pfragmentList;
    }

    public final void addFragment(final Fragment fragment) {
        this.fragmentList.add(fragment);
    }

    public final List<String> getAtomsList() {
        return atomList;
    }

    public final void setAtomList(final List<String> patomList) {
        this.atomList = patomList;
    }

    public final int getAtomListSize() {
        return this.atomList.size();
    }

    public final void addAtom(final String record) {
        this.atomList.add(record);
    }

    public final void removeAtomByIndex(final int idx) {
        this.atomList.remove(idx);
    }
}
