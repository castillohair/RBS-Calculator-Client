"""
Salid RBS Calculator Client

"""

import re

import bs4
import mechanize

# Versions should comply with PEP440.  For a discussion on single-sourcing
# the version across setup.py and the project code, see
# https://packaging.python.org/en/latest/single_source_version.html
__version__ = '0.1.0'

# Mechanize browser
browser = None

def _select_form_by_action(br, action):
    # Go through all forms
    form_index = None
    for i, f in enumerate(br.forms()):
        if f.attrs['action']==action:
            form_index = i
            break
    if form_index is None:
        raise ValueError("form with action {} not found".format(action))
    # Select form
    br.select_form(nr=form_index)

class StartCodon(object):
    def __init__(self, row_tag):
        # Structure of row tag:
        # 0. <td> with empty text
        # 1. <td> with start position
        # 2. <td> with translation initiation rate
        # 3. <td> with Delta G total
        # 4. <td> with Delta G mRNA-rRNA
        # 5. <td> with Delta G spacing
        # 6. <td> with Delta G standby
        # 7. <td> with Delta G start
        # 8. <td> with Delta G mRNA
        # 9. <td> link for mRNA structure
        # 10. <td> with accuracy details and warnings
        tag_children = [c for c in row_tag.find_all('td')]
        self.start_position = int(tag_children[1].text)
        self.tir = float(tag_children[2].text)
        self.DeltaG_total = float(tag_children[3].text)
        self.DeltaG_mRNA_rRNA = float(tag_children[4].text)
        self.DeltaG_spacing = float(tag_children[5].text)
        self.DeltaG_standby = float(tag_children[6].text)
        self.DeltaG_start = float(tag_children[7].text)
        self.DeltaG_mRNA = float(tag_children[8].text)
        
    def __str__(self):
        s = ""
        s += "Start Position: {}\n".format(self.start_position)
        s += "Translation Initiation Rate (au): {}\n".format(self.tir)
        s += "DeltaG_total: {}\n".format(self.DeltaG_total)
        s += "DeltaG_mRNA_rRNA: {}\n".format(self.DeltaG_mRNA_rRNA)
        s += "DeltaG_spacing: {}\n".format(self.DeltaG_spacing)
        s += "DeltaG_standby: {}\n".format(self.DeltaG_standby)
        s += "DeltaG_start: {}\n".format(self.DeltaG_start)
        s += "DeltaG_mRNA: {}".format(self.DeltaG_mRNA)
        return s

class ReverseRBSResult(object):
    def __init__(self, summary_tag, details_tag):
        # Structure of summary tag: one <tr> in a <tbody> with the following:
        # 0. <td> with expand/collapse symbol
        # 1. <td> String with "RBS Rev. Eng."
        # 2. <td> name of job
        # 3. <td> Summary of result, e.g. "8 start codons from 3.50 to 4077.50
        #    T.I.R"
        # 4. <td> Delete button

        # Note: When the job has not been completed, the summary is "Not
        # Completed Yet".
        
        # Structure of details tag: one <tr> in a <tbody> with the following:
        # 0. <tr> with 3 <td>s with the submitted date, CPU time, and organism,
        #    respectively.
        # 1. <tr> with 1 <td> with the algorithm version
        # 2. <tr> with 1 <td> with the string "mRNA sequence"
        # 3. <tr> with the mRNA sequence
        # 4. <tr> with table headings
        # 5-?. one <tr> for each translation initiation site, with the id
        #      "^start_codon_\d+_\d+$"
        # -1. <tr> with the string "All Gibbs free energies..."

        # Note: When the job has not been completed, only the first two rows are
        # present.
        
        # First, verify that this is a reverse calculation
        summary_tag_children = [c for c in summary_tag.find_all('td')]
        if summary_tag_children[1].text != 'RBS Rev. Eng.':
            raise ValueError('tag information does not correspond to reverse '
                'calculation')
        
        # Extract relevant data from summary tag
        self.name = summary_tag_children[2].text
        self.summary = summary_tag_children[3].text
        self.completed = (self.summary != "Not Completed Yet")
        
        # Extract relevant data from details tag
        details_tag_children = [c for c in details_tag.find_all('tr')]
        row1_tags = [c for c in details_tag_children[0].find_all('td')]
        # Remove the leading text "Submitted: "
        self.submitted_date = row1_tags[0].text.strip()[11:]
        # Remove the leading text "CPU Time: "
        self.cpu_time = row1_tags[1].text.strip()[10:]
        # Remove the leading text "Organism: "
        self.organism = row1_tags[2].text.strip()[10:]
        # Remove the leading text "Version: "
        self.version = details_tag_children[1].text.strip()[9:]
        # The following info is only available if the calculation has been
        # completed
        if self.completed:
            # Just remove leading and trailing spaces
            self.sequence = details_tag_children[3].text.strip()
            # Process using start codon class
            self.start_codons = [StartCodon(r) for r in details_tag.find_all(
                id=re.compile('^start_codon_\d+_\d+$'))]
        else:
            self.sequence = None
            self.start_codons = []
        
    def __str__(self):
        s = ""
        if self.completed:
            s += "Name: {}\n".format(self.name)
            s += "Summary: {}\n".format(self.summary)
            s += "Submitted: {}\n".format(self.submitted_date)
            s += "CPU Time: {}\n".format(self.cpu_time)
            s += "Organism: {}\n".format(self.organism)
            s += "Version: {}\n".format(self.version)
            s += "mRNA Sequence: {}\n".format(self.sequence)
            s += "Start codons: (list of {})".format(len(self.start_codons))
        else:
            s += "Name: {}\n".format(self.name)
            s += "(Job has not been completed)\n"
            s += "CPU Time: {}\n".format(self.cpu_time)
            s += "Organism: {}\n".format(self.organism)
            s += "Version: {}".format(self.version)
        return s

class ForwardRBSResult(object):
    def __init__(self, summary_tag, details_tag):
        raise NotImplementedError
        
    def __str__(self):
        raise NotImplementedError

def login(user, password):
    """
    """
    global browser
    browser = mechanize.Browser()
    browser.open("https://salislab.net/software/reverse")
    _select_form_by_action(browser, 'doLogin')
    browser.form['uname'] = user
    browser.form['pwname'] = password
    browser.submit()

def submit_reverse_job(title, mRNA, organism, algorithm_version='v2.0'):
    """
    """
    global browser
    browser.open("https://salislab.net/software/reverse")
    _select_form_by_action(browser, 'doReverseRBS')
    browser.form['title'] = title
    browser.form['mRNA'] = mRNA
    browser.form['organism'] = organism
    browser.form['algorithm_version'] = [algorithm_version,]
    browser.submit()

def get_results():
    """
    """
    global browser
    # Read results website
    results_url = "https://salislab.net/software/Results?method=all&resultsPerPage=50"
    results_html = browser.open(results_url).read()
    # Parse html
    soup = bs4.BeautifulSoup(results_html, 'lxml')
    results = soup.find(id='ResultsContainer')
    # Extract summaries and details
    result_summaries = results.find_all(id=re.compile('^result_summary_\d+$'))
    result_details = results.find_all(id=re.compile('^result_details_\d+$'))
    assert len(result_summaries)==len(result_details)
    # Process results
    results = []
    for s, d in zip(result_summaries, result_details):
        try:
            r = ReverseRBSResult(s, d)
        except ValueError as e:
            r = ForwardRBSResult(s, d)
        results.append(r)

    return results
