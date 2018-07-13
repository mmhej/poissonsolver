
if(typeof(menu) == "string")
{
	document.write('\
			<ul>\
	')


	document.write('<li><a href="' + refdir + 'periodic.html" ')
	if(menu == 'periodic')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">Periodic BCs</a> </li>')



	document.write('\
			</ul>\
	')
}

