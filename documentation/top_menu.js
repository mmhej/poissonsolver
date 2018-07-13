

document.write('\
		<div style="float: left">\
			<table align="center">\
				<tr>\
					<td style="color: #333333;">\
<h1>GREENFISH</h1>\
					</td>\
				</tr>\
				<tr>\
					<td style="color: #9b9b9b; font-size: 18pt; padding-top: 10px;">\
Notes and Documentation\
					</td>\
				</tr>\
			</table>\
		</div>\
		<div style="float: right">\
		</div>\
		<div style="clear: both;">\
		</div>\
')

if(typeof(menu) == "string")
{
	document.write('\
		<div align="center">\
			<ul>\
	')


	document.write('<li><a href="' + refdir + 'index.html" ')
	if(menu == 'index')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">ABOUT</a> </li>')


	document.write('<li><a href="' + refdir + 'notes.html" ')
	if(menu == 'notes')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">NOTES</a> </li>')


	document.write('<li><a href="' + refdir + 'tutorials.html" ')
	if(menu == 'tutorials')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">TUTORIALS</a> </li>')


	document.write('<li><a href="' + refdir + 'program.html" ')
	if(menu == 'program')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">PROGRAM</a> </li>')



	document.write('\
				</ul>\
			</div>\
	')
}



