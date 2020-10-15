suppressMessages(library(argparse))
suppressMessages(library(stringr))
suppressMessages(library(immunarch))
suppressMessages(library(ggplot2))
suppressMessages(library(Cairo))
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--immdata",
                    type="character",
                    default=NULL,
                    help="the mixcr's data")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

width=16
height=12


immdata=readRDS(args$immdata)
print("### Exploratory analysis")
exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("Status"), .meta = immdata$meta)
print(p1)
ggsave(filename =file.path(args$outdir,"exp_vol.pdf"),device = cairo_pdf,height = height,width=width)
dev.off()

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")

p1 <- vis(exp_len)
p2 <- vis(exp_cnt)
p3 <- vis(exp_vol)
print(p1+p2+p3)
ggsave(filename =file.path(args$outdir,"exp_vis.pdf"),device = cairo_pdf,height = height,width=width)
dev.off()


print("### Clonality")
imm_pr <- repClonality(immdata$data, .method = "clonal.prop")
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_rare <- repClonality(immdata$data, .method = "rare")
imm_hom <- repClonality(immdata$data,
  .method = "homeo",
  .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

vis(imm_top) + vis(imm_top, .by = "Status", .meta = immdata$meta)
ggsave(filename =file.path(args$outdir,"Clonality_top.pdf"),device = cairo_pdf,height = height,width=width)
dev.off()

vis(imm_rare) + vis(imm_rare, .by = "Status", .meta = immdata$meta)
ggsave(filename =file.path(args$outdir,"Clonality_Rare.pdf"),device = cairo_pdf,height = height,width=width)
dev.off()


vis(imm_hom) + vis(imm_hom, .by = c("Status"), .meta = immdata$meta)
ggsave(filename =file.path(args$outdir,"Clonality_Hom.pdf"),device = cairo_pdf,height = height,width=width)
dev.off()




