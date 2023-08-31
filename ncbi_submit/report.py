#!/usr/bin/env python3
"""Some useful classes for analyzing report.xml files
"""

from xml.dom import minidom
import io
import re
from typing import Literal


class ActionReader:
    """Stores/retrieves details about a submission action from NCBI's report*.xml files"""

    def __init__(self,submission_id,submission_status,action) -> None:
        """Create an object to store and easily parse details about an action

        Args:
            `submission_id` (_type_): id of associated submission
            `submission_status` (str): status of this action: "processed-ok","processing","queued","submitted","processed-error"
            `action` (str): the xml associated with the action
        """

        self.submission_id = submission_id
        self.submission_status = submission_status
        self.action = action
        self.target_db = self._getTargetDB()
        self.sample_name = self._getSampleName()
        self.status = self._getStatus()
        self.error_source = self._locateError()
        self.bioproject_id = self._get_bioproject_id() if self.submission_status == "submitted" else None
        self.biosample = self._getBioSample()
        self.sra = self._getSra()

    def __str__(self) -> str:
        """str representation of action"""

        if self.biosample:
            return f"Action ({self.sample_name}: db='{self.target_db}' status='{self.status}' bs_acc='{self.biosample}')"
        else:
            return f"Action ({self.sample_name}: db='{self.target_db}' status='{self.status}')"
    
    def _get_bioproject_id(self):
        """Get the bioproject id associated with this action"""

        meta_node = self.action.getElementsByTagName("Meta")
        if len(meta_node) > 0: return meta_node[0].getAttribute("BioProject")
        # return self.action.getElementsByTagName("Meta")[0].getAttribute("BioProject")


    def _getStatus(self):
        """Returns status of action (from Response tag, which is always present)"""

        return self.action.getElementsByTagName("Response")[0].getAttribute("status")

    def _locateError(self):
        """Returns source of error if status="processed-error" (else None)"""

        if self.status == "processed-error": return self.action.getElementsByTagName("Response")[0].getAttribute("error_source")
        else: return None

    def _getTargetDB(self):
        """Returns target_db of action"""

        return self.action.getAttribute("target_db")

    def _getSampleName(self):
        """Returns name of sample (all caps)"""

        action_id = self.action.getAttribute("action_id").upper()
        return action_id.split(f'{self.submission_id}-')[-1].split(f'_{self.target_db.upper()}')[0]

    def _getBioSample(self):
        """Returns BioSample accession if db is biosample"""

        if self.target_db != "BioSample": return None
        return self.action.getElementsByTagName("Object")[0].getAttribute("accession")

    def _getSra(self):
        """Returns SRA accession if db is SRA"""

        if self.target_db != "SRA": return None
        object_node = self.action.getElementsByTagName("Object")
        if len(object_node)==0:
            return None
        return object_node[0].getAttribute("accession")

    def getIsolate(self,ignore_failure=False):
        """Returns isolate name of sample
        
        Args:
            `ignore_failure` (bool): Whether to ignore cases where the isolate can be found or else raise an exception
        """
        
        iso_node = self.action.getElementsByTagName("Isolate")
        if len(iso_node) == 0:
            if ignore_failure: return None
        iso = iso_node[0].toxml()
        tagless = re.sub('<[^<]+>', "", iso)
        return tagless

class Report:
    """Stores/retrieves information about a submission report and it's actions"""

    def __init__(self,report_file) -> None:
        """Create an object to store and parse details about a report*.xml file 

        Args:
            `report_file` (str | Path): path to a report*.xml file
        """

        self.report_file = report_file
        self.dom = minidom.parse(io.open(report_file))
        self.submission_status,self.submission_id = self._getSubmissionDetails()
        self.status_dicts = {}
        # self._createStatusDict()
        self.sra_actions,self.biosample_actions,self.genbank_actions = {},{},{}
        self.system_errors = 0 # number of sra samples where error_source=="system"
        self._addActions()
        self.accession_dict = {}

    def _getSubmissionDetails(self):
        """Returns Tuple(overall_status, submission_id)"""

        submission = self.dom.getElementsByTagName("SubmissionStatus")[0]
        return submission.getAttribute("status"),submission.getAttribute("submission_id")

    def _prep_status_dict_if_needed(self,action):
        """Creates initial dict with empty lists for target_db, if not yet made

        Args:
            `action` (Action): an object holding an individual action's details
        """

        if not action.target_db in self.status_dicts.keys():
            self.status_dicts[action.target_db] = {"processed-ok":[],"processing":[],"queued":[],"submitted":[],"processed-error":[],"deleted":[]}

    def _store_status(self,action):
        """Fills status_dict[target_db] with dict of possible status values if db present in report

        Args:
            `action` (Action): the action object to collect status info on
        """

        self._prep_status_dict_if_needed(action)
        if action.status == "processed-error" and action.target_db=="SRA": self.system_errors += 1
        # add sample_name to list within appropriate dict
        self.status_dicts[action.target_db][action.status].append(action.sample_name)
        
    def _addActions(self):
        """Stores an action object (and status details) for each action in the XML report"""

        action_dicts = {"SRA":self.sra_actions,"BioSample":self.biosample_actions,"GenBank":self.genbank_actions}
        actions = self.dom.getElementsByTagName("Action")
        for xml_action in actions:
            if xml_action.getElementsByTagName("Response"):
                action = ActionReader(submission_id=self.submission_id,submission_status=self.submission_status,action=xml_action)
                action_dicts[action.target_db][action.sample_name] = action
                self._store_status(action)

    def simpleReport(self):
        """Returns status of submission as a whole ("processed-ok","processing","queued","submitted","processed-error")"""

        return self.submission_status

    def statusReport(self,test_dir=False):
        """Returns listed report on numbers of samples of each status and lists failed samples

        Args:
            `test_dir` (bool): Whether this report is for a test submission. Defaults to False.
        """

        report = []
        for db,status_dict in self.status_dicts.items():
            report.append(f"{db} status report")
            report.append("status\t\tsamples")
            num_actions = 0
            for status,samples in status_dict.items():
                num_actions += len(samples)
                report.append(f"{status}\t{len(samples)}")
            num_processed_error = len(status_dict["processed-error"])
            if num_processed_error > 0 and "BioSample" in db:
                report.append("For help with BioSample errors, see:")
                report.append("  https://www.ncbi.nlm.nih.gov/projects/biosample/docs/submission/validation/errors.xml")
            if num_processed_error > 0:
                # if only some samples failed, list them
                if num_actions != num_processed_error:
                    failed:list = status_dict["processed-error"]
                    report.append(f"failed samples: {failed}")
                # if all SRA samples have both `status="processed-error"` and `error_source="system"` --> probably complete
                elif self.system_errors == num_processed_error and test_dir and db=="SRA":
                    report.append("")
                    report.append('All samples have both `status="processed-error"` and `error_source="system"`.')
                    report.append(' This can occur when test SRA submissions are:')
                    report.append('    "processed as well as the test submission')
                    report.append('     database can process Nanopore submissions"')
                    report.append(' and may be indicative of a successful submission.')
        return report

    def submissionComplete(self,db) -> bool:
        """Returns True if (completed samples exist & SRA submissions failed or are still processing), else False
        
        Args:
            `db` (str): db to check. Options: [SRA,BioSample,(TODO: GenBank)]
        """

        if not db in self.status_dicts.keys():
            return False
        bs_ok = len(self.status_dicts[db]["processed-ok"])
        bs_fails = len(self.status_dicts[db]["processed-error"])
        bs_processing = len(self.status_dicts[db]["processing"])
        bs_submitted = len(self.status_dicts[db]["submitted"])
        return sum((bs_fails,bs_processing,bs_submitted)) == 0 and bs_ok > 0

    # def getBiosampleAccessionDict(self,by_sample_name=False,ignore_failure=False) -> dict:
    #     """Returns a dict of sample -> biosample_accession

    #     Args:
    #         `by_sample_name` (bool, optional): A flag to use sample_name as dict key. Defaults to False.
    #             ``True``: key: 'sample_name'
    #             ``False``: key: 'isolate'
    #         `ignore_failure` (bool, optional): A flag to ignore failed samples. Defaults to False.

    #     Raises:
    #         Exception: warns if trying to create accession dict when submission is incomplete
    #     """

    #     if not ignore_failure and not self.submissionComplete("BioSample"):
    #         raise Exception("all BioSamples must have been successfully submitted to get accession_dict. To get accessions anyway, use the `ignore_failure` flag.")
    #     if by_sample_name == False:
    #         self.accession_dict = {action.getIsolate():action.biosample for action in self.biosample_actions.values()}
    #     else:
    #         self.accession_dict = {sample_name:action.biosample for sample_name,action in self.biosample_actions.items()}
    #     return self.accession_dict

    # def getSraAccessionDict(self,by_sample_name=False,ignore_failure=False) -> dict:
    #     """Returns a dict of sample -> sra_accession

    #     Args:
    #         `by_sample_name` (bool, optional): A flag to use sample_name as dict key. Defaults to False.
    #             ``True``: key: 'sample_name'
    #             ``False``: key: 'isolate'
    #         `ignore_failure` (bool, optional): A flag to ignore failed samples. Defaults to False.

    #     Raises:
    #         Exception: warns if trying to create accession dict when submission is incomplete
    #     """

    #     if not ignore_failure and not self.submissionComplete("SRA"):
    #         raise Exception("all SRA samples must have been successfully submitted to get accession_dict. To get accessions anyway, use the `ignore_failure` flag.")
    #     if by_sample_name == False:
    #         self.accession_dict = {action.getIsolate():action.sra for action in self.sra_actions.values()}
    #     else:
    #         self.accession_dict = {sample_name:action.sra for sample_name,action in self.sra_actions.items()}
    #     return self.accession_dict

    def getAccessionDict(self,db:Literal["SRA","BioSample"],by_sample_name=False,ignore_failure=False) -> dict:
        """Returns a dict of sample -> sra_accession

        Args:
            `db` (str): db to check. Options: [SRA,BioSample]
            `by_sample_name` (bool, optional): A flag to use sample_name as dict key. Defaults to False.
                ``True``: key: 'sample_name'
                ``False``: key: 'isolate'
            `ignore_failure` (bool, optional): A flag to ignore failed samples. Defaults to False.

        Raises:
            Exception: warns if trying to create accession dict when submission is incomplete
        """
        # TODO: add GenBank?
        if not self.submissionComplete(db) and not ignore_failure:
            raise Exception(f"all {db} samples must have been successfully submitted to get accession_dict. To get accessions anyway, use the `ignore_failure` flag.")
        db_lower = db.lower()
        if db_lower == "sra": # can't get isolate from here
            by_sample_name = True
        if by_sample_name == False:
            self.accession_dict = {action.getIsolate(ignore_failure=ignore_failure):getattr(action,db_lower) for action in getattr(self,f"{db_lower}_actions").values()}
        else:
            self.accession_dict = {sample_name:getattr(action,db_lower,None) for sample_name,action in getattr(self,f"{db_lower}_actions").items()}
        return self.accession_dict
