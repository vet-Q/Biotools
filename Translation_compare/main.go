package main

import (
	"fmt"
	"strings"
)

var NC string = "mgcrhqlwvg ycqygvhqdl cqklgarkal senkileify nvqyvktssn iilchellsd nplflnnaql klrifgeldt lsinftldni sfnemltryw ysmailyklt eaiqyfyqry shfkdwrlic gvaynnvfdl heiynkektn ididemmqla cmydcnytti yycfmlgadi nramitsvmn fcegnlflci dlgadafees meiasqtnnw ilinillfkn yspdssllsi kttdpekina lldeekyksk nmliyeeslf hiygvni"

var LAO2 string = `mgcrhqlwvg ycqygvhqdl cqklgarkal senkileify nvqyvktssn iilchellsd nplflnnaql klrifgeldt lsinftldni sfnemltryw ysmailyklt eaiqyfyqry dldbfbwj shfkdwrlic gvaynnvfdl heiynkektn ididemmklt cstydgnyst iyycfmlgad inramitsvm nfcegnlflc idlgadafee smeiasqtnn wilinillfk nycpdsslls ikttdpekin alldeekyks knmliyeesl fhiygvni`

func main() {
	NC = strings.ReplaceAll(NC, " ", "")
	LAO2 = strings.ReplaceAll(LAO2, " ", "")
	for i, _ := range NC {
		if NC[i] == LAO2[i] {
			fmt.Printf("OK!, %s,%s\n", string(NC[i]), string(LAO2[i]))
		} else {
			fmt.Printf("NOPE!, %s,%s\n", string(NC[i]), string(LAO2[i]))
		}
	}
}
