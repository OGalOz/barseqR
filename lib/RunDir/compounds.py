# python3
# Python Translation of feba/lib/Compounds.pm
import os
import re
import logging
from RunDir.feba_utils import read_table, read_column_names

# LoadCompounds happens before LoadMedia

# First run
def LoadCompounds(compounds_dir, compounds_dict, synonyms_dict):
    """
    We take compounds info from the file 'compounds.tsv', located in this 
        directory/metadata/Compounds.tsv
        The file has the following headers:
            Compound, CAS, Company, CatalogNumber, FW, Synonyms, Hans80Anti, 
                Hans80metals, FEBA_carbon, FEBA_nitrogen, FEBA_stress, All_star

    Args:
        compounds_dir: (str) Path to dir ('metadata')
        compounds_dict: Starts as empty dict, becomes
            compounds: dict
                compound_name -> [compound_name, CAS (str), MW (str)]
        synonyms_dict: Starts as empty dict, becomes
            synonyms: dict
                synonym_name (str) -> compound_name (str)

    """
    # Uses subfunction SynToKey, read_table, read_column_names 

    compounds_file = os.path.join(compounds_dir, "Compounds.tsv")
    req = "Compound CAS FW Synonyms".split(" ")
    compounds = read_table(compounds_file, req)
    headers = read_column_names(compounds_file)
    if not len(compounds) > 0:
        raise Exception("No rows in " + compounds_file)

    for row in compounds:
        compound = row["Compound"].strip()
        if compound in compounds_dict:
            raise Exception("Duplicate compound id {}".format(compound))
        # mw means molar mass/ molecular weight. CAS is an ID.
        mw = row["FW"]
        if mw == "NA":
            mw = ""
        if not (mw == "" or re.search(r"^[0-9]+[.]?[0-9]*$", mw)):
            raise Exception(
                "Invalid weight {} for compound {} in {}\n \
                        ".format(
                    mw, compound, compounds_file
                )
            )
        cas = row["CAS"].strip()
        if cas == "NA":
            cas = ""

        # Below accounts for cases where cas is multiple CAS divided by '/'
        cas = cas.split("/")[0].strip()
        if not (cas == "" or (re.search(r"^\d+[0-9-]*\d$", cas))):
            raise Exception(
                "Invalid cas number '{}' for compound {} in {}\n \
                        ".format(
                    mw, compound, compounds_file
                )
            )
        compounds_dict[compound] = [compound, cas, mw]

        # Synonyms
        syns = [compound]
        if re.search(r";", row["Synonyms"]):
            # allow semicolon separators instead of comma separators
            syns += row["Synonyms"].split("; ")
        else:
            syns += row["Synonyms"].split(", ")
        for syn in syns:
            key = SynToKey(syn)
            if key == "":
                continue
            if key in synonyms_dict and synonyms_dict[key] != compound:
                logging.warning(
                    "Warning: non-unique synonym "
                    "{}: {} or {}\n".format(syn, compound, synonyms_dict[key])
                )
                synonyms_dict[key] = compound

    return compounds_dict, synonyms_dict 


def SynToKey(syn):
    """ We remove everything BUT the values a-z...+-, and make lower case
        so for example "  Org- acid." becomes "Org-acid" at first,
        then it becomes "org-acid"

    """
    syn = re.sub(r"[^a-zA-Z0-9+-]", '', syn)
    return syn.lower()


# Second Run
# Returns two hashes:
# media => list of components
# media => attributes
# Handles both the media file and the mixes file
# Does NOT look up compound synonyms or otherwise check the results
def LoadMedia(metadir, compounds, synonyms, unknownComponents, reuseComponents):
    """
    Args:
        metadir (dir path)
        compounds: dict
            compound_name -> [compound_name, CAS (str), MW (str)]
        synonyms: dict
            synonym_name (str) -> compound_name (str)
        unknownComponents: dict
            compound_name -> 1 if Unknown
        reuseComponents: dict
            compound -> media_d
                media_d:
                    media_name -> 1 IF some condition met
    Returns:
        media_dict: (d)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mix_dict: (d) (Like media_dict)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mixAttr: (d)
            attribute (str) -> value (str) e.g.
                Description -> Defined minimal media for soil and groundwater bacteria with glucose and MES buffer
                or
                Minimal -> TRUE
        compounds: (d)
            compound_name -> [compound_name, CAS (str), MW (str)]
        synonyms: (d)
            synonym_name (str) -> compound_name (str)
    """
    # Uses subfunctions ParseMediaFile, ExpandMedia, SetupComponentList,
    # SetupExclude
    mediaFile = os.path.join(metadir, "media")
    mixFile = os.path.join(metadir, "mixes")
    if not os.path.isfile(mediaFile):
        raise Exception("No such file: " + mediaFile + "\n")
    if not os.path.isfile(mixFile):
        raise Exception("No such file: " + mixFile + "\n")

    # Init mediaExclude - media => excluded compound => 1
    # This is first occurence of mediaExclude
    mediaExclude = {}
    # The following are dicts
    # media/mix is key => list<list of len(3)>
    media, mediaAttr, mediaExclude = ParseMediaFile(mediaFile, mediaExclude)
    mix, mixAttr, mediaExclude = ParseMediaFile(mixFile, mediaExclude)

    # Validation:
    # Each media must have a Description
    # media should not also be mixes
    for k in mediaAttr.keys():
        med, attr = k, mediaAttr[k]
        if not ("Description" in attr and attr["Description"] != ""):
            raise Exception("No description for media " + med)
        if med in mix:
            raise Exception("Media {} is also a mix".format(med))

    # Each mix must have a Description and a numeric X value
    # mixes should not also be media
    for k in mixAttr.keys():
        mx, attr = k, mixAttr[k]
        if not ("Description" in attr and attr["Description"] != ""):
            raise Exception("No description for mix " + mx)
        if not ("X" in attr and (re.search(r"^[0-9]+[.]?[0-9]*$", attr["X"]))):
            raise Exception("Invalid X for mix " + mx)
        if mx in media:
            raise Exception("Mix {} is also a media".format(mx))

    # Expand media in terms of other media
    nMaxCycles = 100
    for nCycle in range(nMaxCycles + 2):
        if nCycle >= nMaxCycles:
            raise Exception(
                "Too many cycles of media expansion -- is there " "a cycle?\n"
            )
        nChanges = 0
        for md in media.keys():
            nChanges += ExpandMedia(md, media, mediaExclude)[0]
        if nChanges == 0:
            break

    # Replace compound synonyms with compounds, and record any that are not
    # known or are duplicates
    # Below media[k] is a list of length 3
    for k in media.keys():
        media, unknownComponents, reuseComponents = SetupComponentList(
            k,
            media[k],
            media,
            mix,
            unknownComponents,
            reuseComponents,
            compounds,
            synonyms,
        )
    for k in mix.keys():
        media, unknownComponents, reuseComponents = SetupComponentList(
            k,
            mix[k],
            media,
            mix,
            unknownComponents,
            reuseComponents,
            compounds,
            synonyms,
        )

    # Replace compound synonyms with compounds
    # Warn if excluded compound is in the media (and remove it)
    # Define "minus" mixes as needed
    for k in mediaExclude.keys():
        media, mix, mixAttr, unknownComponents = SetupExclude(
            k,
            mediaExclude[k],
            unknownComponents,
            media,
            mix,
            mixAttr,
            compounds,
            synonyms,
        )
    results_dict = {
            "media_dict": media,
            "mix_dict": mix,
            "mixAttr": mixAttr,
            "compounds_dict": compounds,
            "synonyms_dict": synonyms
            }
    return results_dict


def SetupExclude(media, excludeHash, unknownComponents, 
                media_dict, mix, mixAttr, compounds, synonyms):
    """
    Args:
        media: (str) Media name

        excludeHash: (d) value from hashing above media string in dict mediaExclude
                    compound_name (str) -> 1 IF compound included in media

        unknownComponents: (d) 
            compound_name -> 1 if Unknown
        media_dict: (d)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mix: (d) (Like media_dict)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mixAttr: (d)
            attribute (str) -> value (str) e.g.
                Description -> Defined minimal media for soil and groundwater bacteria with glucose and MES buffer
                or
                Minimal -> TRUE
        compounds: (d)
            compound_name -> [compound_name, CAS (str), MW (str)]
        synonyms: (d)
            synonym_name (str) -> compound_name (str)

    """

    # Uses subfunctions FindCompound

    # We set up the list of excluded
    excluded = []
    for orig in excludeHash.keys():
        compound = FindCompound(orig, compounds, synonyms)
        if compound is not None:
            excluded.append(compound)
        else:
            unknownComponents[orig] = 1

    # Dict of excluded compound => number of alterations
    excluded_dict = {x: 0 for x in excluded}

    # And update the media and any incorporated mixes to actually exclude
    updatedComponents = []
    for component in media_dict[media]:
        compound, number, units = component
        if compound in excluded_dict:
            excluded_dict[compound] += 1
        elif compound in mix:
            # compound is a mix -- check if it needs to be altered
            mixExcluded = {}
            mixKeep = []
            for mixComponent in mix[compound]:
                mixCompound = mixComponent[0]
                if mixCompound in excluded_dict:
                    mixExcluded[mixCompound] = 1
                else:
                    mixKeep.append(mixComponent)
            if len(mixExcluded.keys()) > 0:
                # Define the new mix and use it instead
                minusString = " ".join(
                    ["minus " + x for x in sorted(mixExcluded.keys())]
                )
                mixNew = compound + " " + minusString
                if mixNew in mix:
                    if not (len(mix[mixNew])) == (len(mixKeep)):
                        raise Exception(
                            "Wrong number of components for mix "
                            + mixNew
                            + " which can also be built via exclude"
                            " from " + compound
                        )
                else:
                    mix[mixNew] = mixKeep
                    mixAttr[mixNew] = mixAttr[compound]
                    mixAttr[mixNew]["Description"] += " " + minusString
                updatedComponents.append([mixNew, number, units])
                for mixCompound in mixExcluded.keys():
                    excluded_dict[mixCompound] += 1
            else:
                updatedComponents.append(component)
        else:
            updatedComponents.append(component)
    media_dict[media] = updatedComponents
    for k in excluded_dict.keys():
        compound, count = k, excluded_dict[k]
        if count is None:
            raise Exception("B112 " + media + compound)
        if not count > 0:
            logging.warning(
                "Warning: excluding {} from media".format(compound)
                + "{} had no effect\n".format(media)
            )

    return [media_dict, mix, mixAttr, unknownComponents]


def ParseMediaFile(mediaFile, mediaExclude):
    """
    Args:
        mediaFile: str file path
        mediaExclude: (either empty dict or dict as listed below)

    Returns:
        comp: (d)
            media (str) -> list<compound_l>
            where compound_l list<compound (str), concentration (str), units (str)>
            e.g. [Ammonium chloride, 0.25, g/L]
        attr: (d)
            attribute (str) -> value (str) e.g.
                Description -> Defined minimal media for soil and groundwater bacteria with glucose and MES buffer
                or
                Minimal -> TRUE
        mediaExclude: (d)
            media (str) -> excluded compounds (d)
                excluded compounds (d)
                    compound_name (str) -> 1 IF compound included in media

    """
    # Uses subfunction ListValidUnits
    comp = {} # returns as "media or mix" (if media file or mix file)
    attr = {} # returns as "mediaAttr or mixAttr" (like above)
    COMPOUND, NUMBER, UNITS = 0, 1, 2
    validUnits = {x: 1 for x in ListValidUnits()}
    validAttr = {x: 1 for x in ["Description", "Minimal", "X"]}
    curMedia = None
    # readingCompounds is a boolean
    readingCompounds = 0
    with open(mediaFile, "r") as f:
        mf_lines = f.read().split("\n")
    for j in range(len(mf_lines)):
        l = mf_lines[j]
        # Handle DOS mode files
        l = re.sub(r'[\r\n]+$', '', l)
        # strip trailing fields that are empty
        l = re.sub(r'\t+$', '', l)
        F = l.split("\t")
        if len(F) == 1:
            curMedia = None
        elif re.search(r"^#", F[0]):
            pass
        elif len(F) == 2:
            att, val = F
            if att == "Media":
                curMedia = val
                curMedia = re.sub(r" +$", "", curMedia)
                if curMedia in comp:
                    raise Exception("Duplicate media entry for " + curMedia)
                comp[curMedia] = []
                attr[curMedia] = {}
                readingCompounds = 0
            elif att in validAttr:
                if curMedia is None:
                    raise Exception(
                        "No media id yet at line:\n" "{}\n in {}".format(l, mediaFile)
                    )
                if att in attr[curMedia]:
                    raise Exception(
                        "Duplicate attr {} for media {}".format(att, curMedia)
                    )
                attr[curMedia][att] = val
            else:
                raise Exception("Invalid media attribute {}".format(F[0]))
        elif len(F) == 3:
            if curMedia is None:
                raise Exception(
                    "No media id yet at line:\n" "{}\nin {}".format(
                        l, mediaFile)
                    )
            if (
                (re.search(r"^Controlled", F[0]))
                and F[1] == "Concentration"
                and F[2] == "Units"
            ):
                readingCompounds = 1
            else:
                if not readingCompounds == 1:
                    raise Exception(
                        "No compounds header for " "{} at \n{}\n".format(curMedia, l)
                    )
                compound, concentration, units = F
                compound = re.sub(r" +$", "", compound)
                concentration = re.sub(r"^ +", "", concentration)
                concentration = re.sub(r" +$", "", concentration)
                if concentration == "-" and units == "-":
                    if curMedia in mediaExclude:
                        mediaExclude[curMedia][compound] = 1
                    else:
                        mediaExclude[curMedia] = {compound: 1}
                else:
                    cr = concentration
                    if not (
                        cr == ""
                        or re.search(r"^\d+$", cr)
                        or (re.search(r"^\d+[.]\d*$", cr))
                        or (re.search(r"^\d+[.]?\d*[eE][+-]\d+$", cr))
                    ):
                        raise Exception(
                            "Invalid concentration "
                            + concentration
                            + " in line\n"
                            + l
                            + "\nfor "
                            + curMedia
                            + "\n"
                            "in " + mediaFile
                        )
                    units.strip()
                    if units == "":
                        raise Exception("unit missing at line {} in media file: {} ".format(
                                        j, mediaFile))
                    if units not in validUnits:
                        raise Exception(
                            "Invalid unit "
                            + units
                            + " for\n"
                            + l
                            + "\nin "
                            + curMedia
                            + ", "
                            + mediaFile
                        )
                    comp[curMedia].append([compound, cr, units])
        else:
            raise Exception("Wrong number of fields in line:\n"
                    "{}\nIn file {} at line {} ".format(l, mediaFile, j) \
                            + " Num Fields: {}".format(len(F)))
    return [comp, attr, mediaExclude]


def ExpandMedia(media_str, media_dict, mediaExclude):
    """

    Args:
        media_str: (str) Name of the media
        media_dict: (d)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mediaExclude: (d)
            media (str) -> excluded compounds (d)
                excluded compounds (d)
                    compound_name (str) -> 1 IF compound included in media
    We go through the medias and if a media component is
        itself another media, then we expand the media to
        include the subcomponents of that component.
        e.g. Media A contains 10mL of Media B, and Media B
        is composed of compounds d, e, f. So we add d,e,f
        in the right concentrations to Media A
    """
    # Uses no sub functions
    components = []
    nExpand = 0
    for row in media_dict[media_str]:
        # comp must mean component
        comp, conc, units = row
        conc = float(conc)
        if comp == media_str:
            raise Exception(
                "Invalid definition of medium "
                "{} -- it includes itself\n".format(media_str)
            )
        if comp in media_dict:
            # comp is another media
            if comp in mediaExclude:
                raise Exception(
                    "Invalid definition of medium "
                    "{} - it includes medium {} with an exclude".format(media_str, comp)
                    + "compound"
                )
            if not units == "X":
                raise Exception(
                    "Invalid definition of medium "
                    "{} - it includes medium {} but units are not X".format(
                        media_str, comp
                    )
                )
            for subComponent in media_dict[comp]:
                comp2, conc2, units2 = subComponent
                conc2 = float(conc2)
                conc2 *= conc
                components.append([comp2, conc2, units2])
            nExpand += 1
        else:
            components.append(row)
    media_dict[media_str] = components
    return [nExpand, mediaExclude]


# media is a str, l is a list of lists of length 3, everything else dict? 
def SetupComponentList(
    media,
    l,
    media_dict,
    mix_dict,
    unknownComponents,
    reuseComponents,
    compounds,
    synonyms,
):
    """ We take media name and list of compounds within media and
    Args:
        media: (str)
        l: list<compound_l> hashed key media in media_dict like below. The media's compound list
        media_dict: (d)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mix_dict: (d)
            like above media_dict
        unknownComponents: (d) 
            compound_name -> 1 if Unknown
        reuseComponents: dict 
            compound -> media_d
                media_d:
                    media_name -> 1 IF some condition met
        compounds: dict
            compound_name (str) -> [compound_name, CAS (str), MW (str)]
        synonyms: dict
            synonym (str) -> compound_name (str)

    """
    # Uses subfunctions FindCompound
    # 1 or 0
    isMedia = True if media in media_dict else False
    isMix = True if media in mix_dict else False
    if not (isMedia or isMix):
        raise Exception("Unknown media:\n" + media)

    # index of compound within l
    COMPOUND = 0
    # transfer synonyms or record that it is unknown
    for row in l:
        orig, undef, units = row
        origInMix = 1 if orig in mix_dict else 0
        if isMedia == 1 and origInMix == 1:
            if not units == "X":
                raise Exception("Mix must be included with X units. Media: {}".format(
                    media
                    ))
        else:
            compound = FindCompound(orig, compounds, synonyms)
            if compound is not None:
                row[COMPOUND] = compound
            else:
                unknownComponents[row[COMPOUND]] = 1

    # Updating media dict if not updated (?) regarding synonyms.
    media_dict[media] = l

    # Record repeat entries
    seen = {}
    for row in l:
        compound = row[COMPOUND]
        if compound in seen:
            if compound in reuseComponents:
                reuseComponents[compound][media] = 1
            else:
                reuseComponents[compound] = {media:1}
        seen[compound] = 1

    return [media_dict , unknownComponents, reuseComponents]


# Complete - syn is str, compounds and synonyms are dicts 
# Run 4th - after GetMediaComponents
def FindCompound(syn, compounds, synonyms):
    """
    syn: str
    compounds: (d)
        compound_name -> [compound_name, CAS (str), MW (str)]
    synonyms: (d)
        synonym_name (str) -> compound_name (str)
    """
    # Uses subfunction SynToKey
    if syn in compounds:
        return syn
    key = SynToKey(syn)
    if key in synonyms:
        return synonyms[key]
    return None

# Run 3rd
def GetMediaComponents(media_str, media_dict, mixAttr, mix):
    """
    Args:
        media_str: (str) Media name

        media_dict: (d)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]

        mixAttr: (d)
            attribute (str) -> value (str) e.g.
                Description -> Defined minimal media for soil and groundwater bacteria with glucose and MES buffer
                or
                Minimal -> TRUE

        mix: (d)
            like above media_dict

    Returns:
        returns undef for unknown media; otherwise, a  
        list<[compound,concentration,units,mix]>
        where mix is empty (not undefined) unless the compound was 
        included indirectly via a mix
    """
    if media_str not in media_dict:
        return None
    out = []
    for row in media_dict[media_str]:
        # comp, units are str. conc float/int
        comp, conc, units = row
        conc = float(conc)
        if comp in mix:
            if not units == "X":
                raise Exception("Units for mix {} in media {} are not X".format(
                    comp, media_str))
            if not "X" in mixAttr[comp]:
                raise Exception("No X value for mix " + comp)
            rel = conc/float(mixAttr[comp]["X"])
            if not rel > 0 and rel < 1e4:
                raise Exception("Invalid relative X {} for {} in {}".format(
                    rel, comp, media_str))
            for row2 in mix[comp]:
                comp2, conc2, units2 = row2
                conc2 = float(conc2)
                out.append([comp2, conc2*rel, units2, comp])
        else:
            out.append([comp, conc, units, ""])

    return out

# Run last (5th)
def GetMixComponents(mix, mix_dict):
    """
    Args:
        mix: (str)
        mix_dict: (mix - > mix_info like listed in other parts of file.)
    """
    if mix not in mix_dict:
        return None
    else:
        return mix_dict[mix]




def ListValidUnits():
    return "g/L mg/L ug/L M mM uM vol% ml/L X".split(" ")
