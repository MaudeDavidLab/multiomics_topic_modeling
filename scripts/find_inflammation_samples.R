library(dplyr)

all_metadata <- readRDS("Lab/M3/data/all_meta_data.rds")

dim(all_metadata)

vars <- c("Functional.bowel.finding..Biospecimen.", "Stool.frequency..Biospecimen.",
          "GI.issues.this.week..M3...Biospecimen.", "GI.issues.two.weeks.ago..M3...Biospecimen.", "GI.issues.one.month.ago..M3...Biospecimen.",
          "Diarrhea..Biospecimen.", "Constipation..Biospecimen.")

tmp <- all_metadata[, c("Biospecimen.Name", "Host.Name", vars)]
issues <- tmp[tmp$GI.issues.this.week..M3...Biospecimen. & tmp$GI.issues.two.weeks.ago..M3...Biospecimen., ]

View(issues[order(issues$Host.Name), ])
#hosts 086_A , 086_N , 088_A . THe A samples are likely inflammed with diarrhea and GI issues, but 086_N is interesting because it they have chronic constipation while their sibling has chronic diarrhea

weird <- all_metadata[all_metadata$Host.Name %in% c("086_A", "086_N"), ]
View(weird[order(weird$Host.Name), ])
View(weird[ , c("Restaurant.prepared.meals..consumption.frequency...Biospecimen.",
                "Meals.prepared.at.home..consumption.frequency...Biospecimen.",
                "Ready.to.eat.meals..consumption.frequency...Biospecimen.",
                "Host.Name")])

weird$Restaurant.prepared.meals..consumption.frequency...Biospecimen.
weird$Meals.prepared.at.home..consumption.frequency...Biospecimen.
weird$Ready.to.eat.meals..consumption.frequency...Biospecimen.



#find particularly healthy samples?
no_issues <- tmp[!tmp$GI.issues.this.week..M3...Biospecimen. &
                   !tmp$GI.issues.two.weeks.ago..M3...Biospecimen.&
                   !tmp$GI.issues.one.month.ago..M3...Biospecimen.&
                   tmp$Functional.bowel.finding..Biospecimen. == "Tends to have normal bowel function" &
                   tmp$Stool.frequency..Biospecimen. == 1, ]
View(no_issues)
