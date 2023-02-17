Submit-block ::= {
  contact {
    contact {
      name name {
        last "Researcher",
        first "Agood",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "University of Education",
        div "Bioinformatics and Genomics",
        city "Somewhere",
        sub "NC",
        country "USA",
        street "123 Any Street",
        email "aresearch@uni.edu",
        phone "123-456-7890",
        postal-code "24680"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Example",
            first "Seymore",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        },
        {
          name name {
            last "Researcher",
            first "Agood",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "University of Education",
        div "Bioinformatics and Genomics",
        city "Somewhere",
        sub "NC",
        country "USA",
        street "123 Any Street",
        postal-code "24680"
      }
    }
  },
  subtype new
}
Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        "PRJNA553747"
      }
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
