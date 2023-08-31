#!/usr/bin/env python3
"""Classes for converting NCBI TSVs to XML files for submission."""

import datetime
from abc import ABC, abstractmethod
from textwrap import indent,dedent
from typing import Literal
import pandas as pd
from ncbi_submit.helpers import series2dict,get_bioproject_spuid
# from ncbi_submit.ncbi import NCBI # only use for testing (with typing hints) - otherwise causes circular import

class Action(ABC):
    """Converts details about an individual sample to an XML action"""

    def __init__(self,ncbi,indent_by=2) -> None:
        """An object that can generate the text for an xml submission <Action> block
        
        Attributes:
            `ncbi` (NCBI): An object storing information about the submission data, goals, and useful methods
            `indent_by` (int): Number of spaces by which to indent these lines. Defaults to 2.
        """

        ncbi:NCBI
        self.ncbi = ncbi
        self.indent_by = indent_by

    @abstractmethod
    def generate_action_lines(self):
        """Generates individual XML lines to form a single action"""
        pass

class GenBank_Action(Action):
    """Creates XML action to submit all samples from plate to GenBank"""

    def generate_action_lines(self):
        """Generates individual XML lines to form a single action"""

        num = f"_{self.ncbi.attemtpt_num}" if self.ncbi.vary_spuid else ""
        yield indent(dedent(f"""\
            <Action>
              <AddFiles target_db="GenBank">
                <File file_path="genbank.zip">
                  <DataType>genbank-submission-package</DataType>
                </File>
                <Attribute name="wizard">BankIt_SARSCoV2_api</Attribute>
                <Attribute name="auto_remove_failed_seqs">no</Attribute>
                <Identifier>
                  <SPUID spuid_namespace="{self.ncbi.centerAbbr}">{self.ncbi.plate}.sarscov2{num}</SPUID>
                </Identifier>
              </AddFiles>
            </Action>"""),
            self.indent_by*" "
        )

class SRA_BioSample_Action(Action):

    def __init__(self,action_type:Literal["sra","biosample"],ncbi,biosample_link,indent_by=2,attributes={}) -> None:
        """Instantiates SRA_BioSample_Action
        
        Attributes
            `action_type` (str): Database that this action is for.
            `ncbi` (NCBI): An object storing information about the submission data, goals, and useful methods
            `biosample_link` (str): xml tag line linking to BioSample (spuid or accession)
            `indent_by` (int): Number of spaces by which to indent these lines. Defaults to 2.
            `attributes` (dict): Attributes to go in the xml attributes sections
        """

        ncbi:NCBI
        self.action_type = action_type
        self.ncbi = ncbi
        self.biosample_link = biosample_link
        self.indent_by = indent_by
        self.attributes = attributes

    @abstractmethod
    def generate_action_lines(self):
        """Generates individual XML lines (or groups of lines) to form a single action"""
        pass

    def generate_xml_attributes(self):
        """Yields xml text for each attribute in the style of the requested db"""

        for attr_name,attr_value in self.attributes.items():
            # skip empty attributes
            if str(attr_value) == "nan": continue
            elif not str(attr_value).strip(): continue
            # this one goes somewhere special
            if "bioproject_accession" in attr_name: continue
            # drop the time from the date (i.e. " 00:00:00")
            if attr_name == "collection_date": attr_value = str(attr_value).split(" ")[0]
            # add attributes with different spacing for different submissions
            if self.action_type == "sra":
                # these attributes are excluded (since they go elsewhere)
                if "filename" in attr_name or "filetype" in attr_name: continue
                name_of_attr = "name"
            elif self.action_type == "biosample":
                name_of_attr = "attribute_name"
            else: raise NotImplementedError(f"db {self.action_type} not accounted for")
            yield f'<Attribute {name_of_attr}="{attr_name}">{attr_value}</Attribute>'

    def get_bioproject_link(self,accession=None):
        """Returns link to bioproject
        
        Args:
            `accession` (str,optional): The bioproject accession to link the submission to.
        """

        if self.ncbi.bioproject_presets["create_new"]:
            # link to BioProject SPUID
            return get_bioproject_spuid(self.ncbi)
        else:
            # use provided accession or else default from `config_file`
            accession = self.ncbi.get_bioproject_accession(accession)
            # link to known BioProject
            return f'<PrimaryId db="BioProject">{accession}</PrimaryId>'

    def reference_bioproject(self): # this could be replaced by a function that only links biosample if not creating it as a seperate Action step
        """Returns xml text for referencing BioProject accession"""

        return dedent(f"""\
        <AttributeRefId name="BioProject">
          <RefId>
            {self.get_bioproject_link()}
          </RefId>
        </AttributeRefId>""")

class SRA_Action(SRA_BioSample_Action):
    """A class to produce SRA or BioSample XML actions"""

    def __init__(self,ncbi,biosample_link,sra_link,indent_by=2,attributes={},accession_dict=None) -> None:
        """Instantiates SRA_Action
        
        Attributes:
            `ncbi` (NCBI): An object storing information about the submission data, goals, and useful methods
            `biosample_link` (str): xml tag line linking to BioSample (spuid or accession)
            `sra_link` (str): xml tag line containing SRA spuid
            `indent_by` (int): Number of spaces by which to indent these lines. Defaults to 2.
            `attributes` (dict): Attributes to go in the xml attributes sections
            `accession_dict` (dict): maps sample_name -> biosample accession. If sample already has accession, xml will link to accession otherwise, uses spuid
        """

        super().__init__("sra",ncbi=ncbi,biosample_link=biosample_link,indent_by=indent_by,attributes=attributes)
        self.sra_link = sra_link
        self.accession_dict = accession_dict
        self.ncbi:NCBI

    def link_bioSample(self,biosample_spuid):
        """Returns SRA xml text referencing corresponding BioSample

        Args:
          `biosample_spuid` (str): spuid or biosample accession tag line
        """

        return dedent(f"""\
        <AttributeRefId name="BioSample">
          <RefId>
            {biosample_spuid}
          </RefId>
        </AttributeRefId>""")

    def generate_file_info(self):
        """Yields xml text for adding fastq files to an SRA sample"""

        for key,val in self.attributes.items():
            if key.startswith("filename"):
                if str(val) == "nan":
                    raise Exception("Filename looks like 'nan' - verify that all samples are labeled correctly in this plate's gisaid_uploader.log file")
                yield dedent(f"""\
                    <File file_path="{val}">
                      <DataType>generic-data</DataType>
                    </File>""")
            
    def generate_action_lines(self):
        """Generates individual XML lines (or groups of lines) to form a single action"""

        yield indent(dedent(f"""\
            <Action>
              <AddFiles target_db="SRA">"""),
            self.indent_by*" "
        )
        for file_xml in self.generate_file_info():
            yield indent(dedent(file_xml),(4+self.indent_by)*" ")
        
        for attribute in self.generate_xml_attributes():
            yield indent(attribute,"      ")
        # self.generate_xml_attributes(row,db="sra",blankspace="      ").lstrip()
        yield indent(f"""{self.reference_bioproject()}""",(4+self.indent_by)*" ")
        yield indent(self.link_bioSample(self.biosample_link),(4+self.indent_by)*" ")
        yield indent(dedent(f"""\
                <Identifier>
                  {self.sra_link}
                </Identifier>
              </AddFiles>
            </Action>"""),
            self.indent_by*" "
        )

class BioSample_Action(SRA_BioSample_Action):
    """A class to produce SRA XML actions"""

    def __init__(self,ncbi,biosample_link,indent_by=2,attributes={}) -> None:
        """Instantiates BioSample_Action
        
        Attributes
            `ncbi` (NCBI): An object storing information about the submission data, goals, and useful methods
            `biosample_link` (str): xml tag line linking to BioSample (spuid or accession)
            `indent_by` (int): Number of spaces by which to indent these lines. Defaults to 2.
            `attributes` (dict): Attributes to go in the xml attributes sections
        """

        super().__init__("biosample",ncbi=ncbi,biosample_link=biosample_link,indent_by=indent_by,attributes=attributes)

    def generate_action_lines(self):
        """Generates individual XML lines (or groups of lines) for BioSample action"""

        # package options and details: https://www.ncbi.nlm.nih.gov/biosample/docs/packages/
        # SARS-CoV-2.cl.1.0 is assumed, currently but this could be moved to config
        package = self.ncbi.biosample_presets["package"]

        yield indent(dedent(f"""\
            <Action>
              <AddData target_db="BioSample">
                <Data content_type="XML">
                  <XmlContent>
                    <BioSample schema_version="2.0">
                      <SampleId>
                        {self.biosample_link}
                      </SampleId>
                      <Descriptor>
                        <Title>{self.ncbi.plate} BioSample</Title>
                      </Descriptor>
                      <Organism>
                        <OrganismName>{self.ncbi.biosample_presets["organism"]}</OrganismName>
                      </Organism>"""),
            self.indent_by*" "
        )
        row_acc = self.attributes.get("bioproject_accession")
        yield indent(dedent(f"""\
                      <BioProject>
                        <PrimaryId db="BioProject">{self.ncbi.get_bioproject_accession(row_acc)}</PrimaryId>
                      </BioProject>"""),
            (10+self.indent_by)*" "
        )
        yield indent(dedent(f"""\
                      <Package>{package}</Package>
                      <Attributes>"""),
            (10+self.indent_by)*" "
        )
        for attribute in self.generate_xml_attributes():
            yield indent(attribute,"              ")

        yield indent(dedent(f"""\
                      </Attributes>
                    </BioSample>
                  </XmlContent>
                </Data>"""),
            (4+self.indent_by)*" "
        )
        yield indent(dedent(f"""\
                <Identifier>
                  {self.biosample_link}
                </Identifier>
              </AddData>
            </Action>"""),
            self.indent_by*" "
        )        

class Submission(ABC):
    """Converts TSV to XML submission"""

    # @abstractmethod
    def __init__(self,ncbi,data_type:Literal["BS+SRA","BS","SRA","GB"],**kwargs) -> None:
        """Abstract class for creating submission xml files
        
        Attributes:
            `ncbi` (NCBI): An object storing information about the submission data, goals, and useful methods
            `data_type` (str): name of submission type
        """

        ncbi:NCBI
        self.ncbi = ncbi
        self.data_type = data_type
        # for key,val in kwargs.items():
        #     setattr(self,key,val)

    @abstractmethod
    def xml_actions(self):
        """Generates indvidual XML actions or lines of XML actions"""
        pass

    def get_title(self):
        """Returns unique title line for plate"""

        if self.ncbi.plate:
            return f"""<Title>{self.ncbi.plate} - {self.ncbi.biosample_presets["isolation_source"].lower()} - sars-cov-2</Title>"""
        else: return ""

    def xml_description(self):
        """Returns lab-specific header for any NCBI xml file"""

        ## TODO: allow date specifications in `config_file`?
        # next_week = datetime.datetime.now() + datetime.timedelta(days=7)
        # release = next_week.strftime('%Y-%m-%d')
        today = datetime.datetime.now()
        release = today.strftime('%Y-%m-%d') #TODO: add this to `config_file`/parameters
        test = "test " if self.ncbi.test_dir==True else ""
        return indent(dedent(f"""\
            <Description>
              {self.get_title()}
              <Comment>SARS-CoV-2 {test}submission - {self.ncbi.plate} {self.data_type}</Comment>
              <Organization type="center" role="owner">
                <Name>{self.ncbi.affiliation['div']}</Name>
                <Contact email="{self.ncbi.email}">
                  <Name>
                    <First>{self.ncbi.contact[1]}</First>
                    <Last>{self.ncbi.contact[0]}</Last>
                  </Name>
                </Contact>
              </Organization>
              <Hold release_date="{release}"/>
            </Description>"""
        ),"  ")
    
    def generate_xml_lines(self):
        """Yields each piece of XML document"""

        yield '<?xml version="1.0"?>'
        yield "<Submission>"
        yield self.xml_description()
        for line in self.xml_actions(): yield line
        yield "</Submission>"

    def write_xml(self,filename):
        """Combines all XML pieces and writes to '`filename`'"""

        with filename.open("w") as out:
            for line in self.generate_xml_lines():
                end = "" if line.endswith("\n") else "\n"
                out.write(line+end)

class GenBank_Submission(Submission):
    """Creates GenBank submission XML file"""

    def __init__(self,ncbi) -> None:
        super().__init__(ncbi,data_type="GB")

    def xml_actions(self):
        """Generates indvidual XML actions or lines of XML actions"""
        
        for line in GenBank_Action(self.ncbi).generate_action_lines():
            yield line

class SRA_BioSample_Submission(Submission):
    """Converts SRA and BioSample TSVs to submission XML"""

    def __init__(self,ncbi,add_biosample=True,add_sra=True,update_reads=False,spuid_endings=None,report_files=[],download_reports=True) -> None:
        """Instantiate a Submission object
        
        Attributes:
            `ncbi` (NCBI): An object storing information about the submission data, goals, and useful methods
            `add_biosample` (bool): whether to add biosample xml actions to the submission
            `add_sra` (bool): whether to add sra xml actions to the submission
            `update_reads` (bool): whether this submission includes any samples for which the reads are being updated (were previously submitted)
            `spuid_endings` (str | dict): SPUID suffixes to use for specified samples.
              (as str with ":" or ","): 'suffix1:samp1,samp2;suffix2:samp3'. Strings will be converted to the dict format.
              (as str without ":" or ","): all samples will get the str as their suffix
              (as dict): {samp1:suffix1,samp2:suffix1,samp3:suffix2}
            `biosample_df` (DataFrame): df containing final data for biosample_df submission
            `sra_df` (DataFrame): df containing final data for sra_df submission
            `report_files` (list, optional): List of report files from which to seek accessions, if provided. Defaults to None.
            `download_reports` (bool, optional): If True and `update_reads` is True, reports will be downloaded. Defaults to False.
        Auto-set attributes:
            `data_type` (str): name of submission type ("BS+SRA") # TODO: remove?
            `db` (str): database for which to collect report info ("bs_sra")
        """

        ncbi:NCBI
        if update_reads:
            if not ncbi.ftp:
                ncbi.ncbiConnect()
            ncbi._get_all_accessions("bs_sra",report_files=report_files,download_reports=download_reports)
        super().__init__(ncbi,data_type="BS+SRA")
        self.add_biosample = add_biosample
        self.add_sra = add_sra
        self.db = "bs_sra"
        self.update_reads = update_reads
        self.spuid_endings = self.parseSpuidEndings(spuid_endings)

    def parseSpuidEndings(self,spuid_endings=None):
        """Converts user-provided `spuid_endings` to a dictionary

        Args:
            `spuid_endings` (str | dict): SPUID suffixes to use for specified samples.
              (as str with ":" or ","): 'suffix1:samp1,samp2;suffix2:samp3'. Strings will be converted to the dict format.
              (as str without ":" or ","): all samples will get the str as their suffix
              (as dict): {samp1:suffix1,samp2:suffix1,samp3:suffix2}
        """

        if not spuid_endings:
            return {}
        if type(spuid_endings) is dict:
            return spuid_endings
        elif type(spuid_endings) is str:
            spuid_endings = spuid_endings.strip("'\"")
            if not ":" in spuid_endings:
                return spuid_endings
            endings_dict = {}
            suffix_groups = spuid_endings.split(";")
            for group in suffix_groups:
                suffix = group.split(":")[0]
                samples = group.split(":")[1].split(",")
                for sample in samples:
                    endings_dict[sample] = suffix
            return endings_dict
        else: raise TypeError("`spuid_endings` must be a string or a dictionary")

    def getSpuidOrLink(self,row,sra_only,accession_dict):
        """Returns SPUID or link for SRA and BioSample submissions
        
        Args:
            `row` (pd.Series|dict): sample-specific details
            `sra_only` (bool): whether samples are only being submitted to SRA this time
            `accession_dict` (dict): maps sample_name -> biosample accession. If sample already has accession, xml will link to accession otherwise, uses spuid
        """

        num = f"_{self.ncbi.attemtpt_num}" if self.ncbi.vary_spuid else ""
        # print(row["sample_name"])
        suffix = self.spuid_endings.get(row["sample_name"],"") if type(self.spuid_endings)==dict else self.spuid_endings
        sra_link = f"""<SPUID spuid_namespace="{self.ncbi.centerAbbr}">{row["sample_name"]}_SRA{num}{suffix}</SPUID>"""
        if sra_only or (self.update_reads and row["sample_name"] in accession_dict.keys()):
            # assuming BioSamples were previously submitted, link to their accessions, not their name
            biosample_link = f"""<PrimaryId db="BioSample">{accession_dict[row["sample_name"]]}</PrimaryId>"""
            # xml_text = self.sra_action(row[sra_cols],sra_link,biosample_link,bioproject_accession) # temp
        else:
            biosample_link = f"""<SPUID spuid_namespace="{self.ncbi.centerAbbr}">{row["sample_name"]}_BioSample{num}{suffix}</SPUID>"""
            # xml_text = self.biosample_action(row[biosample_cols],biosample_link,plate,test_dir) + "\n" + self.sra_action(row[sra_cols],sra_link,biosample_link,bioproject_accession)
        return sra_link,biosample_link

    def xml_actions(self):
        """Generates indvidual XML actions or portions of XML actions"""

        sra_cols = self.ncbi.sra["final_cols"]
        biosample_cols = self.ncbi.biosample["final_cols"]

        # TODO: test this with actual `test_dir` submission - then make sure BioProject accession can be retrieved from logs
        # create new BioProject, if requested in `config_file`
        if self.ncbi.bioproject_presets["create_new"]:
            yield self.add_BioProject_xml()

        # required if only submitting SRA (and BioSample has already been submitted)
        # TODO: test with `sra_only`
        if self.add_biosample and not self.add_sra:
            accession_dict = self.ncbi.get_accessions(as_dict=True)
        elif self.update_reads:
            acc_df = self.ncbi.sra_bs_accessions_df[["sample_name","BioSample"]]
            acc_df = acc_df[acc_df["sample_name"].isin(self.ncbi.sra["df"]["sample_name"])]
            accession_dict = dict(zip(acc_df["sample_name"],acc_df["BioSample"]))
            if len(accession_dict)==0: raise AttributeError(f"No previously-submitted BioSample accessions found but the flag --update_reads was used. \nTo specify reports containing those accessions, use the --report_files flag. \nAccessions can also be written directly in a csv and pointed to via the variable `extra_accessions` (currently {self.ncbi.extra_accessions}) in your config file (currently {self.ncbi.config_file})")
        else: accession_dict = {}

        # alternate adding BioSample and SRA actions for each sample
        for i,row in self.ncbi.merged["df"].iterrows():
            sra_link,biosample_link = self.getSpuidOrLink(row=row, sra_only=not self.add_biosample, accession_dict=accession_dict)
            if self.add_biosample and not (self.update_reads and row["sample_name"] in accession_dict.keys()):
                for line in BioSample_Action(
                    ncbi=self.ncbi,
                    attributes=series2dict(row,biosample_cols),
                    biosample_link=biosample_link,
                    ).generate_action_lines():
                    yield line
            if self.add_sra:
                for line in SRA_Action(
                    ncbi=self.ncbi,
                    attributes=series2dict(row,sra_cols),
                    biosample_link=biosample_link,
                    sra_link=sra_link,
                    accession_dict=accession_dict,
                    ).generate_action_lines():
                    yield line

    def add_BioProject_xml(self):
        """Returns xml action to create a new BioProject
        For more, see [NCBI BioProject Info](https://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/submit/public-docs/bioproject/bioproject.xsd?view=co)
        
        """
        bioproject_spuid = get_bioproject_spuid(self.ncbi)
        return indent(dedent(f"""\
        <Action>
          <AddData target_db="BioProject">
            <Data content_type="xml">
              <XmlContent>
                <Project schema_version="2.0">
                  <ProjectID>
                    {bioproject_spuid}
                  </ProjectID>
                  <Descriptor>
                    <Title>{self.bioproject_presets["spuid"]}</Title>
                      <Description>
                        <p>{self.bioproject_presets["description"]}</p>
                      </Description>
                    <ExternalLink label="Website: {self.bioproject_presets["website_name"]}">
                      <URL>{self.bioproject_presets["url"]}</URL>
                    </ExternalLink>
                    <Relevance>
                      <Medical>Yes</Medical>
                    </Relevance>
                  </Descriptor>
                  <ProjectType>
                    <ProjectTypeSubmission sample_scope="{self.bioproject_presets["scope"]}">
                      <IntendedDataTypeSet>
                        <DataType>{self.bioproject_presets["dataType"]}</DataType>
                      </IntendedDataTypeSet>
                    </ProjectTypeSubmission>
                  </ProjectType>
                </Project>
              </XmlContent>
            </Data>
            <Identifier>
              {bioproject_spuid}
            </Identifier>
          </AddData>
        </Action>"""
        ),"    ")
